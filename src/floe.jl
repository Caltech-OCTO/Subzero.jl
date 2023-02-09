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
    mc_x::Vector{FT}        # x-coordinates for monte carlo integration
    mc_y::Vector{FT}        # y-coordinates for monte carlo integration
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
        xfloe       <Vector{Float}> vector of floe x-coordinates
        yfloe       <Vector{Float}> vector of floe y-coordinates
        rmax        <Int> floe maximum radius
        area        <Int> floe area
        T           <Float> datatype simulation is run in - either Float64 of Float32
Outputs:
        mc_x        <Vector{T}> vector of monte carlo point x-coordinates
        mc_y        <Vector{T}> vector of monte carlo point y-coordinates
        mc_in       <Vector{Bool}> vector of bools representing if each
                                   monte carlo point is within the given polygon
        err         <Float> error associated with given monte carlo points
"""
function generate_mc_points(npoints, xfloe, yfloe, rmax, area, ::Type{T} = Float64) where T
    mc_x = rmax * (2rand(T, Int(npoints)) .- 1)
    mc_y = rmax * (2rand(T, Int(npoints)) .- 1)
    mc_in = inpoly2(hcat(mc_x, mc_y), hcat(xfloe, yfloe))
    err = abs(sum(mc_in)/npoints * 4 * rmax^2 - area)/area
    return mc_x, mc_y, mc_in, err
end

"""
    Floe(poly::LG.Polygon, hmean, Δh, ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, t::Type{T} = Float64)

Constructor for floe with LibGEOS Polygon
Inputs:
        poly        <LibGEOS.Polygon> 
        h_mean      <Real> mean height for floes
        Δh          <Real> variability in height for floes
        grid        <Grid>
        ρi          <Real> ice density kg/m3 - default 920
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
        mc_n        <Real> number of monte carlo points
        t           <Float> datatype to run simulation with - either Float32 or Float64
Output:
        Floe with needed fields defined - all default field values used so all forcings start at 0 and floe is "alive".
        Velocities and the density of ice can be optionally set.
Note:
        Types are specified at Float64 below as type annotations given that when written LibGEOS could exclusivley use Float64 (as of 09/29/22).
        When this is fixed, this annotation will need to be updated.
        We should only run the model with Float64 right now or else we will be converting the Polygon back and forth all of the time. 
"""
function Floe(poly::LG.Polygon, hmean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n = 1000.0, t::Type{T} = Float64) where T
    floe = rmholes(poly)
    centroid = LG.GeoInterface.coordinates(LG.centroid(floe))::Vector{Float64}
    h = hmean + (-1)^rand(0:1) * rand() * Δh  # floe height
    area = LG.area(floe)::Float64  # floe area
    mass = area * h * ρi  # floe mass
    coords = LG.GeoInterface.coordinates(floe)::PolyVec{Float64}
    moment = calc_moment_inertia(coords, centroid, h, ρi = ρi)
    origin_coords = translate(coords, -centroid)
    ox, oy = seperate_xy(origin_coords)
    rmax = sqrt(maximum([sum(c.^2) for c in origin_coords[1]]))
    angles = calc_poly_angles(coords, T)
    # Generate Monte Carlo Points
    count = 1
    alive = true
    mc_x, mc_y, mc_in, err = generate_mc_points(mc_n, ox, oy, rmax, area, T)
    while err > 0.1
        mc_x, mc_y, mc_in, err = generate_mc_points(mc_n, ox, oy, rmax, area, T)
        count += 1
        if count > 10
            err = 0.0
            alive = false
        end
    end
    mc_x = mc_x[mc_in[:, 1] .|  mc_in[:, 2]]
    mc_y = mc_y[mc_in[:, 1] .|  mc_in[:, 2]]

    return Floe(centroid = convert(Vector{T}, centroid), coords = convert(PolyVec{T}, coords),
                height = convert(T, h), area = convert(T, area), mass = convert(T, mass),
                rmax = convert(T, rmax), moment = convert(T, moment), angles = angles,
                u = convert(T, u), v = convert(T, v), ξ = convert(T, ξ),
                mc_x = mc_x, mc_y = mc_y, alive = alive)
end

"""
    Floe(coords::PolyVec{Float64}, h_mean, Δh, ρi = 920.0, u = 0.0,
    v = 0.0, ξ = 0.0, t::Type{T} = Float64) where T

Floe constructor with PolyVec{Float64}(i.e. Vector{Vector{Vector{Float64}}}) coordinates
Inputs:
        coords      <Vector{Vector{Vector{Float64}}}> floe coordinates
        h_mean      <Real> mean height for floes
        Δh          <Real> variability in height for floes
        grid        <Grid>
        rho_ice     <Real> ice density kg/m3
        u           <Real> x-velocity of the floe - default 0.0
        v           <Real> y-velcoity of the floe - default 0.0
        ksi         <Real> angular velocity of the floe - default 0.0
        t           <Type> datatype to convert ocean fields - must be a Float!
Output:
        Floe with needed fields defined - all default field values used so all forcings and velocities start at 0 and floe is "alive"
"""
Floe(coords::PolyVec{<:Real}, h_mean, Δh; ρi = 920.0, u = 0.0, v = 0.0, ξ = 0.0, mc_n = 1000, t::Type{T} = Float64) where T =
    Floe(LG.Polygon(convert(PolyVec{Float64}, valid_polyvec!(rmholes(coords)))), h_mean, Δh; ρi = ρi, u = u, v = v, ξ = ξ, mc_n = mc_n, t = T) 
    # Polygon convert is needed since LibGEOS only takes Float64 - when this is fixed convert can be removed


function initialize_floe_field(coords::Vector{PolyVec}, h_mean, Δh; ρi = 920.0, mc_n = 1000, t::Type{T} = Float64) where T
    floe_arr = StructArray([Floe(c, h_mean, Δh, ρi = ρi, mc_n = mc_n, t = T) for c in coords])
    # Initialize floe IDs
    floe_arr.id .= range(1, length(floe_arr))
    return floe_arr
end


"""
nfloes is a proxy for floe size. It is now many floes you would want to fit in the domain at your given concentrations without any topography.
Topography decreases the number of floes that will be created, so increase n floes accordingly 
"""
function initialize_floe_field(nfloes::Int, concentrations, domain, h_mean, Δh; ρi = 920.0, mc_n = 1000, t::Type{T} = Float64) where T
    floe_arr = StructArray{Floe{T}}(undef, 0)
    # Split domain into cells with given concentrations
    nrows, ncols = size(concentrations[:, :])
    rowlen = (domain.north.val - domain.south.val) / nrows
    collen = (domain.east.val - domain.west.val) / ncols
    # Availible space in whole domain
    open_water = LG.Polygon(rect_coords(domain.west.val, domain.east.val, domain.south.val, domain.north.val))
    if !isempty(domain.topography)
        open_water = LG.difference(open_water, LG.MultiPolygon(domain.topography.coords))
    end
    open_water_area = LG.area(open_water)
    # Loop over cells
    for i in range(1, nrows)
        for j in range(1, ncols)
            c = concentrations[i, j]
            if c > 0
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
                xpoints = rand(T, ncell)
                ypoints = rand(T, ncell)
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
                    if LG.area(floe_poly) > 0
                        regions = LG.getGeometries(floe_poly)
                        while !isempty(regions)
                            r = pop!(regions)
                            if !hashole(r)
                                floe = Floe(r, h_mean, Δh, ρi = ρi, mc_n = mc_n, t = T)
                                push!(floe_arr, floe)
                                floes_area += floe.area
                            else
                                region_bottom, region_top = split_polygon_hole(r, T)
                                append!(regions, region_bottom)
                                append!(regions, region_top)
                            end
                        end
                    end
                end
            end
        end
    end
    # Initialize floe IDs
    floe_arr.id .= range(1, length(floe_arr))
    return floe_arr
end
