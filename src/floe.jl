"""
Structs and functions used to define floes and floe fields within Subzero
"""

"""
Singular sea ice floe with fields describing current state.
"""
@kwdef mutable struct Floe{FT<:AbstractFloat}
    # Physical Properties -------------------------------------------------
    centroid::Vector{FT}    # center of mass of floe (might not be in floe!)
    coords::PolyVec{FT}     # floe coordinates
    height::FT              # floe height (m)
    area::FT                # floe area (m^2)
    mass::FT                # floe mass (kg)
    rmax::FT                # distance of vertix farthest from centroid (m)
    moment::FT              # mass moment of intertia
    angles::Vector{FT}      # interior angles of floe
    # Monte Carlo Points ---------------------------------------------------
    mc_x::Vector{FT}        # x-coordinates for monte carlo integration centered at origin
    mc_y::Vector{FT}        # y-coordinates for monte carlo integration centered at origin
    # Velocity/Orientation -------------------------------------------------
    α::FT = 0.0             # floe rotation from starting position in radians
    u::FT = 0.0             # floe x-velocity
    v::FT = 0.0             # floe y-velocity
    ξ::FT = 0.0             # floe angular velocity
    # Status ---------------------------------------------------------------
    alive::Bool = true      # floe is still active in simulation
    id::Int = 0             # floe id - set to index in floe array at start of sim
    ghost_id::Int = 0       # ghost id - if floe is a ghost, ghost_id > 0 representing which ghost it is
    ghosts::Vector{Int} = Vector{Int}()  # indices of ghost floes of given floe
    # Forces/Collisions ----------------------------------------------------
    fxOA::FT = 0.0          # force from ocean and atmos in x direction
    fyOA::FT = 0.0          # force from ocean and atmos in y direction
    trqOA::FT = 0.0         # torque from ocean and Atmos
    hflx_factor::FT = 0.0   # heat flux factor can be multiplied by floe height to get the heatflux
    overarea::FT = 0.0      # total overlap with other floe
    collision_force::Matrix{FT} = [0.0 0.0] 
    collision_trq::FT = 0.0
    interactions::Matrix{FT} = zeros(0, 7)
    # Previous values for timestepping  -------------------------------------
    p_dxdt::FT = 0.0        # previous timestep x-velocity
    p_dydt::FT = 0.0        # previous timestep y-velocity
    p_dudt::FT = 0.0        # previous timestep x-acceleration
    p_dvdt::FT = 0.0        # previous timestep x-acceleration
    p_dξdt::FT = 0.0        # previous timestep time derivative of ξ
    p_dαdt::FT = 0.0        # previous timestep ξ
end

"""
Enum to index into floe interactions field with more intuituve names
"""
@enum InteractionFields begin
    floeidx = 1
    xforce = 2
    yforce = 3
    xpoint = 4
    ypoint = 5
    torque = 6
    overlap = 7
end
"""
Index into interactions field with InteractionFields enum objects
"""
Base.to_index(s::InteractionFields) = Int(s)

"""
    generate_mc_points(npoints, xfloe, yfloe, rmax, area, ::Type{T} = Float64)

Generate monte carlo points, determine which are within the given floe, and the error associated with the points
Inputs:
        npoints     <Int> number of points to generate
        xfloe       <Vector{Float}> vector of floe x-coordinates centered on the origin
        yfloe       <Vector{Float}> vector of floe y-coordinates centered on the origin
        rmax        <Int> floe maximum radius
        area        <Int> floe area
        alive       <Bool> true if floe is alive (i.e. will continue in the simulation), else false
        rng         <RNG> random number generator to generate monte carlo points
        T           <Float> datatype simulation is run in - either Float64 of Float32
Outputs:
        mc_x        <Vector{T}> vector of monte carlo point x-coordinates that are within floe
        mc_y        <Vector{T}> vector of monte carlo point y-coordinates that are within floe
        alive       <Bool> true if floe is alive (i.e. will continue in the simulation), else false
Note: You will not end up with npoints. This is the number originally generated, but any not in the floe are deleted.
The more oblong the floe shape, the less points. 
"""
function generate_mc_points(npoints, xfloe, yfloe, rmax, area, alive, rng, ::Type{T} = Float64) where T
    count = 1
    err = T(1)
    mc_x = zeros(T, npoints)
    mc_y = zeros(T, npoints)
    mc_in = fill(false, npoints)
    while err > 0.1
        if count > 10
            err = 0.0
            alive = false
        else
            mc_x .= rmax * (2rand(rng, T, Int(npoints)) .- 1)
            mc_y .= rmax * (2rand(rng, T, Int(npoints)) .- 1)
            in_on = inpoly2(hcat(mc_x, mc_y), hcat(xfloe, yfloe))
            mc_in .= in_on[:, 1] .|  in_on[:, 2]
            err = abs(sum(mc_in)/npoints * 4 * rmax^2 - area)/area
            count += 1
        end
    end

    return mc_x[mc_in], mc_y[mc_in], alive
end

"""
    Floe(poly::LG.Polygon, hmean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n = 1000.0, t::Type{T} = Float64)

Constructor for floe with LibGEOS Polygon
Inputs:
        poly        <LibGEOS.Polygon> 
        hmean      <Real> mean height for floes
        Δh          <Real> variability in height for floes
        grid        <Grid> simulation grid
        ρi          <Real> ice density kg/m3 - default 920
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ξ         <Real> angular velocity of the floe - default 0.0
        mc_n        <Real> number of monte carlo points
        rng         <RNG> random number generator to generate random floe attributes -
                          default is RNG using Xoshiro256++ algorithm
        t           <Float> datatype to run simulation with - either Float32 or Float64
Output:
        Floe with needed fields defined - all default field values used so all forcings start at 0 and floe is "alive".
        Velocities and the density of ice can be optionally set.
Note:
        Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22).
        When this is fixed, this annotation will need to be updated.
        We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function Floe(poly::LG.Polygon, hmean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n::Int = 1000, rng = Xoshiro(), t::Type{T} = Float64) where T
    floe = rmholes(poly)
    # Floe physical properties
    centroid = LG.GeoInterface.coordinates(LG.centroid(floe))::Vector{Float64}
    height = hmean + (-1)^rand(rng, 0:1) * rand(rng) * Δh
    area = LG.area(floe)::Float64
    mass = area * height * ρi
    coords = LG.GeoInterface.coordinates(floe)::PolyVec{Float64}
    moment = calc_moment_inertia(coords, centroid, height, ρi = ρi)
    angles = calc_poly_angles(coords, T)
    origin_coords = translate(coords, -centroid)
    rmax = sqrt(maximum([sum(c.^2) for c in origin_coords[1]]))
    alive = true
    # Generate Monte Carlo Points
    ox, oy = seperate_xy(origin_coords)
    mc_x, mc_y, alive = generate_mc_points(mc_n, ox, oy, rmax, area, alive, rng, T)
    

    return Floe(centroid = convert(Vector{T}, centroid), coords = convert(PolyVec{T}, coords),
                height = convert(T, height), area = convert(T, area), mass = convert(T, mass),
                rmax = convert(T, rmax), moment = convert(T, moment), angles = angles,
                u = convert(T, u), v = convert(T, v), ξ = convert(T, ξ),
                mc_x = mc_x, mc_y = mc_y, alive = alive)
end

"""
    Floe(coords::PolyVec, hmean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n = 1000, t::Type{T} = Float64)

Floe constructor with PolyVec{Float64}(i.e. Vector{Vector{Vector{Float64}}}) coordinates
Inputs:
        coords      <Vector{Vector{Vector{Float64}}}> floe coordinates
        hmean      <Real> mean height for floes
        Δh          <Real> variability in height for floes
        grid        <Grid> simulationg grid
        ρi          <Real> ice density kg/m3 - default 920
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ξ           <Real> angular velocity of the floe - default 0.0
        mc_n        <Real> number of monte carlo points
        rng         <RNG> random number generator to generate random floe attributes -
                          default is RNG using Xoshiro256++ algorithm
        t           <Float> datatype to run simulation with - either Float32 or Float64
Output:
        Floe with needed fields defined - all default field values used so all forcings and velocities start at 0 and floe is "alive"
"""
Floe(coords::PolyVec, hmean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n::Int = 1000, rng = Xoshiro(), t::Type{T} = Float64) where T =
    Floe(LG.Polygon(convert(PolyVec{Float64}, valid_polyvec!(rmholes(coords)))), hmean, Δh; ρi = ρi, u = u, v = v, ξ = ξ, mc_n = mc_n, rng = rng, t = T) 
    # Polygon convert is needed since LibGEOS only takes Float64 - when this is fixed convert can be removed

"""
    poly_to_floes(floe_poly, ρi, min_floe_area = 0.0)

Split a given polygon into regions and split around any holes before turning each region with an area
greater than the minimum floe size into a floe.
Inputs:
        floe_poly       <LibGEOS.Polygon or LibGEOS.MultiPolygon> polygon/multipolygon to turn onto floes
        hmean          <Float>             average floe height
        Δh              <Float>             height range - floes will range in height from hmean - Δh to hmean + Δh
        ρi              <Float> ice density
        mc_n            <Int> number of monte carlo points
        rng             <RNG> random number generator to generate random floe attributes -
                              default is RNG using Xoshiro256++ algorithm
        min_floe_area   <Float> minimum area for floe creation - default is 0
        t               <Float> datatype to run simulation with - either Float32 or Float64
Output:
        Vector{Floe} vector of floes making up input polygon(s) with area above given minimum floe area.
        Floe polygons split around holes if needed. 
"""
function poly_to_floes(floe_poly, hmean, Δh, ρi, mc_n::Int, rng, min_floe_area = 0, t::Type{T} = Float64) where T
    floes = StructArray{Floe{T}}(undef, 0)
    regions = LG.getGeometries(floe_poly)
    while !isempty(regions)
        r = pop!(regions)
        if LG.area(r) > min_floe_area
            if !hashole(r)
                floe = Floe(r, hmean, Δh, ρi = ρi, mc_n = mc_n, rng = rng, t = T)
                push!(floes, floe)
            else
                region_bottom, region_top = split_polygon_hole(r, T)
                append!(regions, region_bottom)
                append!(regions, region_top)
            end
        end
    end
    return floes
end

"""
    initialize_floe_field(coords::Vector{PolyVec}, domain, hmean, Δh; min_floe_area::T = -1.0, ρi::T = 920.0, mc_n = 1000, t::Type{T} = Float64)

Create a field of floes from a list of polygon coordiantes. User is wanrned if floe's do not meet minimum size requirment. 
Inputs:
        coords          <Vector{PolyVec}>   list of polygon coordinates to make into floes
        domain          <Domain>            model domain 
        hmean          <Float>             average floe height
        Δh              <Float>             height range - floes will range in height from hmean - Δh to hmean + Δh
        min_floe_area   <Float>             if a floe below this minimum floe size is created program will throw a warning (optional) -
                                            default is 0, but if a negative is provided it will be replaced with 4*Lx*Ly/1e4
                                            where Lx and Ly are the size of the domain edges
        ρi              <Float>             ice density (optional) - default is 920.0
        mc_n            <Int>               number of monte carlo points to intially generate for each floe (optional) - 
                                            default is 1000 - note that this is not the number you will end up with as some will be outside of the floe
        rng             <RNG>               random number generator to generate random floe attributes -
                                            default is RNG using Xoshiro256++ algorithm
        T               <Type>              An abstract float type to run the simulation in (optional) - default is Float64
Output:
        floe_arr <StructArray> list of floes created from given polygon coordinates
"""
function initialize_floe_field(coords::Vector{PolyVec{T}}, domain, hmean, Δh; min_floe_area = 0.0, ρi = 920.0, mc_n::Int = 1000, rng = Xoshiro(), t::Type{T} = Float64) where T
    floe_arr = StructArray{Floe{T}}(undef, 0)
    floe_polys = [LG.Polygon(valid_polyvec!(c)) for c in coords]
    # Remove overlaps with topography
    if !isempty(domain.topography)
        topo_poly = LG.MultiPolygon(domain.topography.coords)
        floe_polys = [LG.difference.(f, topo_poly) for f in floes]
    end
    # Turn polygons into floes
    for p in floe_polys
        append!(floe_arr, poly_to_floes(p, hmean, Δh, ρi, mc_n, rng, min_floe_area, T))
    end
    # Warn about floes with area less than minimum floe size
    min_floe_area = min_floe_area > 0 ? min_floe_area : T(4 * (domain.east.val - domain.west.val) * (domain.north.val - domain.south.val) / 1e4)
    if any(floe_arr.area .< min_floe_area)
        @warn "Some user input floe areas are less than the suggested minimum floe size."
    end
    # Warn about floes with centroids outside of domain
    if any(domain.west.val .< first.(floe_arr.centroid) .< domain.east.val) || any(domain.south.val .< last.(floe_arr.centroid) .< domain.north.val)
        @warn "Some floe centroids are out of the domain."
    end
    # Initialize floe IDs
    floe_arr.id .= range(1, length(floe_arr))
    return floe_arr
end

"""
    initialize_floe_field(nfloes::Int, concentrations, domain, hmean, Δh; min_floe_area = -1.0, ρi = 920.0, mc_n::Int = 1000, t::Type{T} = Float64)

Create a field of floes using Voronoi Tesselation.
Inputs:
        nfloes          <Int>       number of floes to try to create - note you might not end up with this number of floes -
                                    topography in domain and multiple concentrations can decrease number of floes created
        concentrations  <Matrix>    matrix of concentrations to fill domain. If size(concentrations) = N, M then split the
                                    domain into NxM cells, each to be filled with the corresponding concentration. 
                                    If concentration is below 0, it will default to 0. If it is above 1, it will default to 1
        domain          <Domain>    model domain 
        hmean          <Float>     average floe height
        Δh              <Float>     height range - floes will range in height from hmean - Δh to hmean + Δh
        min_floe_area   <Float>     if a floe below this minimum floe size is created it will be deleted (optional) -
                                    default is 0, but if a negative is provided it will be replaced with 4*Lx*Ly/1e4
                                    where Lx and Ly are the size of the domain edges
        ρi              <Float>     ice density (optional) - default is 920.0
        mc_n            <Int>       number of monte carlo points to intially generate for each floe (optional) - 
                                    default is 1000 - note that this is not the number you will end up with as some will be outside of the floe
        rng             <RNG>       random number generator to generate random floe attributes -
                                    default is RNG using Xoshiro256++ algorithm
        T               <Type>      An abstract float type to run the simulation in (optional) - default is Float64
Output:
        floe_arr <StructArray> list of floes created using Voronoi Tesselation of the domain with given concentrations.
"""
function initialize_floe_field(nfloes::Int, concentrations, domain, hmean, Δh; min_floe_area = 0.0, ρi = 920.0, mc_n::Int = 1000, rng = Xoshiro(), t::Type{T} = Float64) where T
    floe_arr = StructArray{Floe{T}}(undef, 0)
    # Split domain into cells with given concentrations
    nrows, ncols = size(concentrations[:, :])
    Lx = domain.east.val - domain.west.val
    Ly = domain.north.val - domain.south.val
    rowlen = Ly / nrows
    collen = Lx / ncols
    # Availible space in whole domain
    open_water = LG.Polygon(rect_coords(domain.west.val, domain.east.val, domain.south.val, domain.north.val))
    if !isempty(domain.topography)
        open_water = LG.difference(open_water, LG.MultiPolygon(domain.topography.coords))
    end
    open_water_area = LG.area(open_water)
    min_floe_area = min_floe_area >= 0 ? min_floe_area : T(4 * Lx * Lx / 1e4)
    # Loop over cells
    for i in range(1, nrows)
        for j in range(1, ncols)
            c = concentrations[i, j]
            if c > 0
                c = c > 1 ? 1 : c
                # Grid cell bounds
                xmin = domain.west.val + collen * (j - 1)
                ymin = domain.south.val + rowlen * (i - 1)
                cell_bounds = rect_coords(xmin, xmin + collen, ymin, ymin + rowlen)
                trans_vec = [xmin, ymin]
                # Open water in cell
                open_cell = LG.intersection(LG.Polygon(cell_bounds), open_water)
                open_coords = LG.GeoInterface.coordinates(open_cell)::PolyVec{T}
                open_area = LG.area(open_cell)::T
                # Create points to seed floes
                ncell = ceil(Int, nfloes * open_area / open_water_area / c)
                xpoints = rand(rng, T, ncell)
                ypoints = rand(rng, T, ncell)
                in_points = inpoly2(hcat(collen * xpoints .+ trans_vec[1], rowlen * ypoints .+ trans_vec[2]), reduce(vcat, open_coords))
                in_idx = in_points[:, 1] .|  in_points[:, 2]
                # Create floes
                tess_floes = voronoicells(xpoints[in_idx], ypoints[in_idx], Rectangle(Point2(0.0, 0.0), Point2(1.0, 1.0))).Cells 
                floes_area = T(0.0)
                floe_idx = shuffle(range(1, length(tess_floes)))
                while !isempty(floe_idx) && floes_area/open_area <= c
                    idx = pop!(floe_idx)
                    floe_coords = [valid_ringvec!([Vector(f) .* [collen, rowlen] .+ trans_vec for f in tess_floes[idx]])]
                    floe_poly = LG.intersection(LG.Polygon(floe_coords), open_cell)
                    floes = poly_to_floes(floe_poly, hmean, Δh, ρi, mc_n, rng, min_floe_area, T)
                    append!(floe_arr, floes)
                    floes_area += sum(floes.area)
                end
            end
        end
    end
    # Initialize floe IDs
    floe_arr.id .= range(1, length(floe_arr))
    return floe_arr
end
