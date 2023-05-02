"""
Structs and functions used to define floes and floe fields within Subzero
"""

"""
Enum to index into floe interactions field with more intuituve names
"""
@enum StatusTag begin
    active = 1
    remove = 2
    fuse = 3
end

mutable struct Status
    tag::StatusTag
    fuse_idx::Vector{Int}
end

Status() = Status(active, Vector{Int}())  # active floe

"""
    StressCircularBuffer{FT<:AbstractFloat}

Extended circular buffer for the stress history that hold 2x2 matrices of stress
values and allows for efficently taking the mean of  buffer by keeping an
element-wise running total of values within circular buffer
"""
mutable struct StressCircularBuffer{FT<:AbstractFloat}
    cb::CircularBuffer{Matrix{FT}}
    total::Matrix{FT}
end
"""
    StressCircularBuffer{FT}(capacity::Int)

Create a stress buffer with given capacity
Inputs:
    capacity    <Int> capacity of circular buffer
Outputs:
    StressCircularBuffer with given capacity and a starting total that is a 2x2
    matrix of zeros.
"""
StressCircularBuffer{FT}(capacity::Int) where {FT} =
    StressCircularBuffer{FT}(
        CircularBuffer{Matrix{FT}}(capacity),
        zeros(FT, 2, 2)
    )
"""
    push!(scb::StressCircularBuffer, data)

Adds element to the back of the circular buffer and overwrite front if full.
Add data to total and remove overwritten value from total if full.
Inputs:
    scb     <StressCircularBuffer> stress circular buffer
    data    <Matrix> 2x2 stress data
Outputs:
    Add data to the buffer and update the total to reflect the addition
"""
function Base.push!(scb::StressCircularBuffer, data)
    if scb.cb.length == scb.cb.capacity
        scb.total .-= scb.cb[1]
    end
    scb.total .+= data
    push!(scb.cb, data)
end

"""
fill!(scb::StressCircularBuffer, data)

Grows the buffer up-to capacity, and fills it entirely. It doesn't overwrite
existing elements. Adds value of added items to total.
Inputs:
    scb     <StressCircularBuffer> stress circular buffer
    data    <Matrix> 2x2 stress data
Outputs:
    Fill all empty buffer slots with data and update total to reflect additions
"""
function Base.fill!(scb::StressCircularBuffer, data)
    scb.total .+= (scb.cb.capacity - scb.cb.length) * data
    fill!(scb.cb, data)
end

"""
    mean(scb::StressCircularBuffer)

Calculates mean of buffer, over the capacity of the buffer. If the buffer is not
full, empty slots are counted as zeros.
Inputs:
    scb     <StressCircularBuffer> stress circular buffer
Outputs:
    mean of stress circular buffer over the capacity of the buffer
"""
function Statistics.mean(scb::StressCircularBuffer)
    return scb.total / capacity(scb.cb) 
end

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
    angles::Vector{FT}      # interior angles of floe in degrees
    # Monte Carlo Points ---------------------------------------------------
    mc_x::Vector{FT}        # x-coordinates for monte carlo integration centered
                            #   at origin
    mc_y::Vector{FT}        # y-coordinates for monte carlo integration centered
                            #   at origin
    # Velocity/Orientation -------------------------------------------------
    α::FT = 0.0             # floe rotation from starting position in radians
    u::FT = 0.0             # floe x-velocity
    v::FT = 0.0             # floe y-velocity
    ξ::FT = 0.0             # floe angular velocity
    # Status ---------------------------------------------------------------
    status::Status = Status() # floe is active in simulation
    id::Int = 0             # floe id - set to index in floe array at start of
                            #   sim - unique to all floes
    ghost_id::Int = 0       # ghost id - if floe is a ghost, ghost_id > 0
                            #   representing which ghost it is
                            #   if floe is not a ghost, ghost_id = 0
    parent_ids::Vector{Int} = Vector{Int}()  # if the floe was originally part
                            # of one or several floes, list parent ids
    ghosts::Vector{Int} = Vector{Int}()  # indices of ghost floes of given floe
    # Forces/Collisions ----------------------------------------------------
    fxOA::FT = 0.0          # force from ocean and atmos in x direction
    fyOA::FT = 0.0          # force from ocean and atmos in y direction
    trqOA::FT = 0.0         # torque from ocean and Atmos
    hflx_factor::FT = 0.0   # heat flux factor can be multiplied by floe height
                            #   to get the heatflux
    overarea::FT = 0.0      # total overlap with other floe
    collision_force::Matrix{FT} = zeros(1, 2)
    collision_trq::FT = 0.0
    interactions::Matrix{FT} = zeros(0, 7)
    stress::Matrix{FT} = zeros(2, 2)
    stress_history::StressCircularBuffer{FT} = StressCircularBuffer(1000)
    strain::Matrix{FT} = zeros(2, 2)
    # Previous values for timestepping  -------------------------------------
    p_dxdt::FT = 0.0        # previous timestep x-velocity
    p_dydt::FT = 0.0        # previous timestep y-velocity
    p_dudt::FT = 0.0        # previous timestep x-acceleration
    p_dvdt::FT = 0.0        # previous timestep y-acceleration
    p_dξdt::FT = 0.0        # previous timestep time derivative of ξ
    p_dαdt::FT = 0.0        # previous timestep angular-velocity
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
Create a range of interactions field columns with InteractionFields enum objects
"""
Base.:(:)(a::InteractionFields, b::InteractionFields) = Int(a):Int(b)

"""
    generate_mc_points(
        npoints,
        xfloe,
        yfloe,
        rmax,
        area,
        status,
        rng,
        ::Type{T} = Float64,
    ) where T

Generate monte carlo points, determine which are within the given floe, and the
error associated with the points
Inputs:
    npoints <Int> number of points to generate
    coords  <PolyVec{AbstractFloat}> PolyVec of floe coords centered on origin
    yfloe   <Vector{Float}> vector of floe y-coordinates centered on the origin
    rmax    <Int> floe maximum radius
    area    <Int> floe area
    status  <Status> floe status (i.e. active, fuse, etc in simulation)
    rng     <RNG> random number generator to generate monte carlo points
Outputs:
    mc_x   <Vector{T}> vector of monte carlo point x-coords that are within floe
    mc_y   <Vector{T}> vector of monte carlo point y-coords that are within floe
    status <Status> tag is `active` if points created correctly, else `remove``
Note:
    You will not end up with npoints. This is the number originally generated,
    but any not in the floe are deleted. The more oblong the floe shape, the
    less points. 
"""
function generate_mc_points(
    npoints::Int,
    coords::PolyVec{T},
    rmax::T,
    area::T,
    status,
    rng,
) where {T<:AbstractFloat}
    count = 1
    err = T(1)
    mc_x = zeros(T, npoints)
    mc_y = zeros(T, npoints)
    mc_in = fill(false, npoints)
    while err > 0.1
        if count > 10
            err = 0.0
            status.tag = remove
        else
            mc_x .= rmax * (2rand(rng, T, Int(npoints)) .- 1)
            mc_y .= rmax * (2rand(rng, T, Int(npoints)) .- 1)
            mc_in = points_in_poly(hcat(mc_x, mc_y), coords)
            err = abs(sum(mc_in)/npoints * 4 * rmax^2 - area)/area
            count += 1
        end
    end

    return mc_x[mc_in], mc_y[mc_in], status
end

"""
    Floe(
        poly::LG.Polygon,
        hmean,
        Δh;
        ρi = 920.0,
        u = 0.0,
        v = 0.0,
        ξ = 0.0,
        mc_n = 1000.0,
        t::Type{T} = Float64,
    )

Constructor for floe with LibGEOS Polygon
Inputs:
    poly  <LibGEOS.Polygon> 
    hmean <Real> mean height for floes
    Δh    <Real> variability in height for floes
    grid  <Grid> simulation grid
    ρi    <Real> ice density kg/m3 - default 920
    u     <Real> x-velocity of the floe - default 0.0
    v     <Real> y-velcoity of the floe - default 0.0
    ξ     <Real> angular velocity of the floe - default 0.0
    mc_n  <Real> number of monte carlo points
    rng   <RNG> random number generator to generate random floe attributes
            default is RNG using Xoshiro256++ algorithm
    t     <Float> datatype to run simulation with - either Float32 or Float64
Output:
    <Floe> with needed fields defined - all default field values used so all
        forcings start at 0 and floe's status is "active". Velocities and the
        density of ice can be optionally set.
Note:
    Types are specified at Float64 below as type annotations given that when
    written LibGEOS could exclusivley use Float64 (as of 09/29/22). When this is
    fixed, this annotation will need to be updated. We should only run the model
    with Float64 right now or else we will be converting the Polygon back and
    forth all of the time. 
"""
function Floe(
    poly::LG.Polygon,
    hmean,
    Δh;
    ρi = 920.0,
    u = 0.0,
    v = 0.0,
    ξ = 0.0,
    mc_n::Int = 1000,
    nhistory = 1000,
    rng = Xoshiro(),
    t::Type{T} = Float64,
) where T
    floe = rmholes(poly)
    # Floe physical properties
    centroid = find_poly_centroid(floe)
    height = hmean + (-1)^rand(rng, 0:1) * rand(rng, T) * Δh
    area_tot = LG.area(floe)::Float64
    mass = area_tot * height * ρi
    coords = find_poly_coords(floe)
    moment = calc_moment_inertia(coords, centroid, height, ρi = ρi)
    angles = calc_poly_angles(coords)
    translate!(coords, -centroid[1], -centroid[2])
    rmax = sqrt(maximum([sum(c.^2) for c in coords[1]]))
    status = Status()
    # Generate Monte Carlo Points
    mc_x, mc_y, status = generate_mc_points(
        mc_n,
        coords,
        rmax,
        area_tot,
        status,
        rng,
    )
    translate!(coords, centroid[1], centroid[2])
    # Generate Stress History
    stress_history = StressCircularBuffer{T}(nhistory)
    fill!(stress_history, zeros(T, 2, 2))

    return Floe(
        centroid = convert(Vector{T}, centroid),
        coords = convert(PolyVec{T}, coords),
        height = convert(T, height),
        area = convert(T, area_tot),
        mass = convert(T, mass),
        rmax = convert(T, rmax),
        moment = convert(T, moment),
        angles = angles,
        u = convert(T, u),
        v = convert(T, v),
        ξ = convert(T, ξ),
        mc_x = mc_x,
        mc_y = mc_y,
        stress_history = stress_history,
        status = status,
    )
end

"""
    Floe(
        coords::PolyVec,
        hmean,
        Δh;
        ρi = 920.0,
        u = 0.0,
        v = 0.0,
        ξ = 0.0,
        mc_n = 1000,
        t::Type{T} = Float64,
    )

Floe constructor with PolyVec{AbstractFloat} coordinates
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
    rng         <RNG> random number generator to generate random floe attributes
                    - default is RNG using Xoshiro256++ algorithm
    t           <Float> datatype to run simulation with - either Float32 or 64
Output:
    <Floe> with needed fields defined - all default field values used so all
        forcings and velocities start at 0 and floe's status is "active"
"""
Floe(
    coords::PolyVec,
    hmean,
    Δh;
    ρi = 920.0,
    u = 0.0,
    v = 0.0,
    ξ = 0.0,
    mc_n::Int = 1000,
    nhistory = 1000,
    rng = Xoshiro(),
    t::Type{T} = Float64,
) where T =
    Floe( # Polygon convert is needed since LibGEOS only takes Float64 - when this is fixed convert can be removed
        LG.Polygon(convert(PolyVec{Float64},
        valid_polyvec!(rmholes(coords)))),
        hmean,
        Δh;
        ρi = ρi,
        u = u,
        v = v,
        ξ = ξ,
        mc_n = mc_n,
        nhistory = nhistory,
        rng = rng,
        t = T,
    ) 

"""
    poly_to_floes(
        floe_poly,
        hmean,
        Δh;
        ρi = 920.0,
        u = 0.0,
        v = 0.0,
        ξ = 0.0,
        mc_n::Int = 1000,
        nhistory::Int = 1000,
        rng = Xoshiro(),
        min_floe_area = 0,
    )

Split a given polygon into regions and split around any holes before turning
each region with an area greater than the minimum floe area into a floe.
Inputs:
    floe_poly       <LibGEOS.Polygon or LibGEOS.MultiPolygon> polygon/
                        multipolygon to turn into floes
    hmean           <AbstratFloat> average floe height
    Δh              <AbstratFloat> height range - floes will range in height
                        from hmean - Δh to hmean + Δh
    ρi              <AbstratFloat> ice density
    u               <AbstratFloat> floe u velocity
    v               <AbstratFloat> floe v velocity
    ξ               <AbstratFloat> floe angular velocity
    mc_n            <Int> number of monte carlo points
    nhistory        <Int> number of element in floe's stress history
    rng             <RNG> random number generator to generate random floe
                        attributes - default is RNG using Xoshiro256++ algorithm
    min_floe_area   <AbstratFloat> minimum area for floe creation - default is 0
Output:
    <StructArray{Floe}> vector of floes making up input polygon(s) with area
        above given minimum floe area. Floe's with holes split around holes. 
"""
function poly_to_floes(
    floe_poly,
    hmean::FT,
    Δh::FT;
    ρi = 920.0,
    u = 0.0,
    v = 0.0,
    ξ = 0.0,
    mc_n::Int = 1000,
    nhistory::Int = 1000,
    rng = Xoshiro(),
    min_floe_area = 0,
) where {FT <: AbstractFloat}
    floes = StructArray{Floe{FT}}(undef, 0)
    regions = LG.getGeometries(floe_poly)::Vector{LG.Polygon}
    while !isempty(regions)
        r = pop!(regions)
        a = LG.area(r)
        if a >= min_floe_area && a > 0
            if !hashole(r)
                floe = Floe(
                    r::LG.Polygon,
                    hmean,
                    Δh,
                    ρi = ρi,
                    u = u,
                    v = v,
                    ξ = ξ,
                    mc_n = mc_n,
                    nhistory = nhistory,
                    rng = rng,
                    t = FT,
                )
                push!(floes, floe)
            else
                region_bottom, region_top = split_polygon_hole(r, FT)
                append!(regions, region_bottom)
                append!(regions, region_top)
            end
        end
    end
    return floes
end

"""
    initialize_floe_field(
        coords::Vector{PolyVec{T}},
        domain,
        hmean,
        Δh;
        min_floe_area = 0.0,
        ρi = 920.0,
        mc_n::Int = 1000,
        nhistory::Int = 1000,
        rng = Xoshiro(),
        t::Type{T} = Float64,
    ) where T

Create a field of floes from a list of polygon coordiantes. User is wanrned if
floe's do not meet minimum size requirment. 
Inputs:
    coords          <Vector{PolyVec}> list of polygon coords to make into floes
    domain          <Domain> model domain 
    hmean           <Float> average floe height
    Δh              <Float> height range - floes will range in height from hmean
                        - Δh to hmean + Δh
    min_floe_area   <Float> if a floe below this minimum floe size is created
                        program will throw a warning (optional) -  default is 0
                        where Lx and Ly are the size of the domain edges
    ρi              <Float> ice density (optional) - default is 920.0
    mc_n            <Int> number of monte carlo points to intially generate for
                        each floe (optional) - default is 1000 - note that this
                        is not the number you will end up with as some will be
                        outside of the floe
    rng             <RNG> random number generator to generate random floe
                        attributes - default is RNG using Xoshiro256++ algorithm
    T               <Type> An abstract float type to run the simulation in
                        (optional) - default is Float64
Output:
    floe_arr <StructArray> list of floes created from given polygon coordinates
"""
function initialize_floe_field(
    coords::Vector{PolyVec{T}},
    domain,
    hmean,
    Δh;
    min_floe_area = 0.0,
    ρi = 920.0,
    mc_n::Int = 1000,
    nhistory::Int = 1000,
    rng = Xoshiro(),
    t::Type{T} = Float64,
) where T
    floe_arr = StructArray{Floe{T}}(undef, 0)
    floe_polys = [LG.Polygon(valid_polyvec!(c)) for c in coords]
    # Remove overlaps with topography
    if !isempty(domain.topography)
        topo_poly = LG.MultiPolygon(domain.topography.coords)
        floe_polys = [LG.difference(f, topo_poly) for f in floe_polys]
    end
    # Turn polygons into floes
    for p in floe_polys
        append!(
            floe_arr, 
            poly_to_floes(
                p,
                T(hmean),
                T(Δh);
                ρi = ρi,
                mc_n = mc_n,
                nhistory = nhistory,
                rng = rng,
                min_floe_area = min_floe_area,
            ),
        )
    end
    # Warn about floes with area less than minimum floe size
    min_floe_area = min_floe_area > 0 ?
        min_floe_area :
        T(
            4 * (domain.east.val - domain.west.val) *
            (domain.north.val - domain.south.val) / 1e4
        )
    if any(floe_arr.area .< min_floe_area)
        @warn "Some user input floe areas are less than the suggested minimum \
            floe area."
    end
    # Warn about floes with centroids outside of domain
    if !all(
        domain.west.val .< first.(floe_arr.centroid) .< domain.east.val) &&
        !all(domain.south.val .< last.(floe_arr.centroid) .< domain.north.val
    )
        @warn "Some floe centroids are out of the domain."
    end
    # Initialize floe IDs
    floe_arr.id .= range(1, length(floe_arr))
    return floe_arr
end

"""
    generate_voronoi_coords(
        desired_points,
        scale_fac,
        trans_vec,
        domain_coords,
        rng;
        max_tries = 10,
        t::Type{T} = Float64,
    ) where T

Generate voronoi coords within a bounding box defined by its lower left corner
and its height and width. Attempt to generate `npieces` cells within the box.
Inputs:
    desired_points  <Int> desired number of voronoi cells
    scale_fac       <Vector{AbstractFloat}> width and height of bounding box -
                        formatted as [w, h] 
    trans_vec       <Vector{AbstractFloat}> lower left corner of bounding box -
                        formatted as [x, y] 
    domain_coords   <Vector{PolyVec{AbstractFloat}}> multipolygon that will
                        eventually be filled with/intersected with the voronoi
                        cells - such as topography
    rng             <RNG> random number generator to generate voronoi cells
    min_to_warn     <Int> minimum number of points to warn if not generated to
                        seed voronoi
    max_tries       <Int> number of tires to generate desired number of points
                        within domain_coords to seed voronoi cell creation
    T               <Type> An abstract float type to run the simulation in
                        (optional) - default is Float64
Outputs:
    coords  <Vector{PolyVec{Float}}> vector of polygon coordinates generated by
        voronoi tesselation. These polygons all fall within the space defined by
        the domain_coords. If less polygons than min_to_warn are generated, the
        user will be warned. 
"""
function generate_voronoi_coords(
    desired_points::Int,
    scale_fac,
    trans_vec,
    domain_coords::Vector{<:PolyVec{<:T}},
    rng,
    min_to_warn::Int;
    max_tries::Int = 10,
    t::Type{T} = Float64,
) where T
    xpoints = Vector{T}()
    ypoints = Vector{T}()
    area_frac = LG.area(LG.MultiPolygon(domain_coords)) / reduce(*, scale_fac)
    # Increase the number of points based on availible percent of bounding box
    npoints = ceil(Int, desired_points / area_frac)
    current_points = 0
    tries = 0
    while current_points < desired_points && tries <= max_tries
        x = rand(rng, T, npoints)
        y = rand(rng, T, npoints)
        # Scaled and translated points
        st_xy = hcat(
            scale_fac[1] * x .+ trans_vec[1],
            scale_fac[2] * y .+ trans_vec[2]
        )
        # Check which of the points are within the domain coords
        in_idx = points_in_poly(st_xy, domain_coords)
        current_points += sum(in_idx)
        tries += 1
        append!(xpoints, x[in_idx])
        append!(ypoints, y[in_idx])
    end
    # If we generated too many cells, remove extra
    if current_points > desired_points
        xpoints = xpoints[1:desired_points]
        ypoints = ypoints[1:desired_points]
        current_points = desired_points
    end
    # Warn if we didn't generate enough cells
    if current_points < min_to_warn
        @warn "Only $current_points floes were able to be generated in \
            $max_tries tries during voronoi tesselation."
    end
    # Make voronoi cells into floes
    coords =
        if current_points > 1
            tess_cells = voronoicells(
                xpoints,
                ypoints,
                Rectangle(Point2(0.0, 0.0), Point2(1.0, 1.0)),
                rng = rng
            ).Cells
            # Scale and translate voronoi coordinates
            [[valid_ringvec!([
                Vector(c) .* scale_fac .+ trans_vec for c in tess
            ])] for tess in tess_cells]
        else
            Vector{Vector{Vector{T}}}()
        end
    return coords
end

"""
    initialize_floe_field(
        nfloes::Int,
        concentrations,
        domain,
        hmean,
        Δh;
        min_floe_area = 0.0,
        ρi = 920.0,
        mc_n::Int = 1000,
        t::Type{T} = Float64)

Create a field of floes using Voronoi Tesselation.
Inputs:
        nfloes          <Int> number of floes to try to create - note you
                            might not end up with this number of floes -
                            topography in domain and multiple concentrations can
                            decrease number of floes created
        concentrations  <Matrix> matrix of concentrations to fill domain. If
                            size(concentrations) = N, M then split the domain
                            into NxM cells, each to be filled with the
                            corresponding concentration. If concentration is
                            below 0, it will default to 0. If it is above 1, it
                            will default to 1
        domain          <Domain> model domain 
        hmean           <Float> average floe height
        Δh              <Float> height range - floes will range in height from
                            hmean - Δh to hmean + Δh
        min_floe_area   <Float> if a floe below this minimum floe size is
                            created it will be deleted (optional) - default is
                            0, but if a negative is provided it will be replaced
                            with 4*Lx*Ly/1e4 where Lx and Ly are the size of the
                            domain edges
        ρi              <Float> ice density (optional) - default is 920.0
        mc_n            <Int> number of monte carlo points to intially generate
                            for each floe (optional) - default is 1000 - note
                            that this is not the number you will end up with as
                            some will be outside of the floe
        rng             <RNG> random number generator to generate random floe
                            attributes - default is RNG using Xoshiro256++
        T               <Type> an abstract float type to run the simulation in
                            (optional) - default is Float64
Output:
        floe_arr <StructArray> list of floes created using Voronoi Tesselation
            of the domain with given concentrations.
"""
function initialize_floe_field(
    nfloes::Int,
    concentrations,
    domain,
    hmean,
    Δh;
    min_floe_area = 0.0,
    ρi = 920.0,
    mc_n::Int = 1000,
    nhistory::Int = 1000,
    rng = Xoshiro(),
    t::Type{T} = Float64
) where T
    floe_arr = StructArray{Floe{T}}(undef, 0)
    # Split domain into cells with given concentrations
    nrows, ncols = size(concentrations[:, :])
    Lx = domain.east.val - domain.west.val
    Ly = domain.north.val - domain.south.val
    rowlen = Ly / nrows
    collen = Lx / ncols
    # Availible space in whole domain
    open_water = LG.Polygon(rect_coords(
        domain.west.val,
        domain.east.val,
        domain.south.val,
        domain.north.val
    ))
    if !isempty(domain.topography)
        open_water = LG.difference(
            open_water, 
            LG.MultiPolygon(domain.topography.coords)
        )
    end
    open_water_area = LG.area(open_water)
    min_floe_area = min_floe_area >= 0 ? min_floe_area : T(4 * Lx * Lx / 1e4)
    # Loop over cells
    for j in range(1, ncols)
        for i in range(1, nrows)
            c = concentrations[i, j]
            if c > 0
                c = c > 1 ? 1 : c
                # Grid cell bounds
                xmin = domain.west.val + collen * (j - 1)
                ymin = domain.south.val + rowlen * (i - 1)
                cell_bounds = rect_coords(
                    xmin,
                    xmin + collen,
                    ymin,
                    ymin + rowlen,
                )
                trans_vec = [xmin, ymin]
                # Open water in cell
                open_cell = LG.intersection(LG.Polygon(cell_bounds), open_water)
                open_coords = find_multipoly_coords(open_cell)
                open_area = LG.area(open_cell)::T
                # Generate coords with voronoi tesselation and make into floes
                ncells = ceil(Int, nfloes * open_area / open_water_area / c)
                floe_coords = generate_voronoi_coords(
                    ncells,
                    [collen, rowlen],
                    trans_vec,
                    open_coords,
                    rng,
                    ncells,
                    t = T,
                )
                nfloes = length(floe_coords)
                if nfloes > 0
                    floes_area = T(0.0)
                    floe_idx = shuffle(rng, range(1, nfloes))
                    while !isempty(floe_idx) && floes_area/open_area <= c
                        idx = pop!(floe_idx)
                        floe_poly = LG.intersection(
                            LG.Polygon(floe_coords[idx]),
                            open_cell
                        )
                        floes = poly_to_floes(
                            floe_poly,
                            T(hmean),
                            T(Δh);
                            ρi = ρi,
                            mc_n = mc_n,
                            nhistory = nhistory,
                            rng = rng,
                            min_floe_area = min_floe_area,
                        )
                        append!(floe_arr, floes)
                        floes_area += sum(floes.area)
                    end
                end
            end
        end
    end
    # Initialize floe IDs
    floe_arr.id .= range(1, length(floe_arr))
    return floe_arr
end
