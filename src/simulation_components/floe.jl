"""
Structs and functions used to define floes and floe fields within Subzero
"""

"""
Enum for differnt floe status
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
    x_subfloe_points::Vector{FT}        # x-coordinates for monte carlo integration centered
                            #   at origin
    y_subfloe_points::Vector{FT}        # y-coordinates for monte carlo integration centered
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
    num_inters::Int = 0
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
    Floe(::Type{FT}, args...; kwargs...)

A float type FT can be provided as the first argument of any Floe constructor. A
Floe of type FT will be created by passing all other arguments to the correct
constructor. 
"""
Floe(::Type{FT}, args...; kwargs...) where {FT <: AbstractFloat} =
    Floe{FT}(args...; kwargs...)

"""
    Floe(args...)

If a type isn't specified, Floe will be of type Float64 and the correct
constructor will be called with all other arguments.
"""
Floe(args...; kwargs...) = Floe{Float64}(args...; kwargs...)

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
    Floe{FT}(
        poly::Polys,
        hmean,
        Δh;
        floe_settings = FloeSettings(),
        rng = Xoshiro(),
        kwargs...
    )

Constructor for floe with LibGEOS Polygon
Inputs:
    poly                <LibGEOS.Polygon> 
    hmean               <Real> mean height for floes
    Δh                  <Real> variability in height for floes
    floe_settings       <FloeSettings> settings needed to initialize floe        
    rng                 <RNG> random number generator to generate floe
                            attributes - default is Xoshiro256++ algorithm
    kwargs      Any other floe fields to set as keyword arguments
Output:
    <Floe> with needed fields defined - all default field values used so all
        forcings start at 0 and floe's status is "active" as long as monte carlo
        points were able to be generated.
"""
function Floe{FT}(
    poly::Polys,
    hmean,
    Δh;
    floe_settings = FloeSettings(),
    rng = Xoshiro(),
    kwargs...
) where {FT <: AbstractFloat}
    floe = rmholes(poly)
    # Floe physical properties
    centroid = collect(GO.centroid(floe))
    height = clamp(
        hmean + (-1)^rand(rng, 0:1) * rand(rng, FT) * Δh,
        floe_settings.min_floe_height,
        floe_settings.max_floe_height,
    )
    area_tot = GO.area(floe)
    mass = area_tot * height * floe_settings.ρi
    coords = find_poly_coords(floe)
    coords = [orient_coords(coords[1])]
    moment = _calc_moment_inertia(FT, floe, centroid, height; ρi = floe_settings.ρi)
    angles = GO.angles(GI.Polygon(coords))
    translate!(coords, -centroid[1], -centroid[2])
    rmax = sqrt(maximum([sum(c.^2) for c in coords[1]]))
    status = Status()
    # Generate Monte Carlo Points
    x_subfloe_points, y_subfloe_points, status = generate_subfloe_points(
        floe_settings.subfloe_point_generator,
        coords,
        rmax,
        area_tot,
        status,
        rng,
    )
    translate!(coords, centroid[1], centroid[2])
    # Generate Stress History
    stress_history = StressCircularBuffer{FT}(floe_settings.nhistory)
    fill!(stress_history, zeros(FT, 2, 2))

    return Floe{FT}(;
        centroid = centroid,
        coords = coords,
        height = height,
        area = area_tot,
        mass = mass,
        rmax = rmax,
        moment = moment,
        angles = angles,
        x_subfloe_points = x_subfloe_points,
        y_subfloe_points = y_subfloe_points,
        stress_history = stress_history,
        status = status,
        kwargs...
    )
end

"""
    Floe{FT}(
        coords::PolyVec,
        hmean,
        Δh;
        ρi = 920.0,
        mc_n = 1000,
        nhistory = 1000,
        rng = Xoshiro(),
        kwargs...,
    )

Floe constructor with PolyVec coordinates
Inputs:
    coords              <Vector{Vector{Vector{Float64}}}> floe coordinates
    hmean               <Real> mean height for floes
    Δh                  <Real> variability in height for floes
    floe_settings       <FloeSettings> settings needed to initialize floe
    rng                 <RNG> random number generator to generate random floe
                            attributes - default uses Xoshiro256++ algorithm
    kwargs              Any other floe fields to set as keyword arguments
Output:
    <Floe> with needed fields defined - all default field values used so all
    forcings start at 0 and floe's status is "active" as long as monte carlo
    points were able to be generated.
"""
Floe{FT}(
    coords::PolyVec,
    hmean,
    Δh;
    floe_settings = FloeSettings(),
    rng = Xoshiro(),
    kwargs...,
) where {FT <: AbstractFloat} =
    Floe{FT}( # Polygon convert is needed since LibGEOS only takes Float64
        LG.Polygon(
            convert(
                PolyVec{Float64},
                valid_polyvec!(rmholes(coords)),
            ),
        ),
        hmean,
        Δh;
        floe_settings = floe_settings,
        rng = rng,
        kwargs...,
    ) 

"""
    poly_to_floes!(
        ::Type{FT},
        floes,
        poly,
        hmean,
        Δh,
        rmax;
        floe_settings,
        rng = Xoshiro(),
        kwargs...
    )

Split a given polygon around any holes before turning each region with an area greater than
the minimum floe area into a floe.
Inputs:
    Type{FT}            <AbstractFloat> Type for grid's numberical fields -
                        determines simulation run type
    floes               <StructArray{Floe}> vector of floes to add new floes to
    poly                <Polygon> polygons to turn into floes
    hmean               <AbstratFloat> average floe height
    Δh                  <AbstratFloat> height range - floes will range in height
                        from hmean - Δh to hmean + Δh
    rmax                <AbstractFloat> maximum radius of floe (could be larger given context)
    floe_settings       <FloeSettings> settings needed to initialize floe
                            settings
    rng                 <RNG> random number generator to generate random floe
                            attributes - default uses Xoshiro256++ algorithm
    kwargs...           Any additional keywords to pass to floe constructor
"""
function poly_to_floes!(
    ::Type{FT},
    floes,
    poly,
    hmean,
    Δh,
    rmax;
    floe_settings = FloeSettings(min_floe_area = 0),
    rng = Xoshiro(),
    kwargs...
) where {FT <: AbstractFloat}
    a = GO.area(poly)
    if a >= floe_settings.min_floe_area && a > 0
        if !hashole(poly)
            floe = Floe(
                FT,
                poly::Polys,
                hmean,
                Δh;
                floe_settings = floe_settings,
                rng = rng,
                kwargs...
            )
            push!(floes, floe)
            return 1
        else
            cx, cy = GO.centroid(GI.gethole(poly, 1), FT)
            new_regions = GO.cut(poly, GI.Line([(cx - rmax, cy), (cx + rmax, cy)]), FT)
            n = 0
            for r in new_regions
                n += poly_to_floes!(FT, floes, r, hmean, Δh, rmax;
                    floe_settings = floe_settings, rng = rng, kwargs...)
            end
            return n
        end
    end
    return 0
end

"""
    initialize_floe_field(args...)

If a type isn't specified, the list of Floes will each be of type Float64 and
the correct constructor will be called with all other arguments.
"""
initialize_floe_field(args...; kwargs...) =
    initialize_floe_field(Float64, args...; kwargs...)

"""
    initialize_floe_field(
        ::Type{FT},
        coords,
        domain,
        hmean,
        Δh;
        floe_settings,
        rng,
    )

Create a field of floes from a list of polygon coordiantes. User is wanrned if
floe's do not meet minimum size requirment. 
Inputs:
    Type{FT}            <AbstractFloat> Type for grid's numberical fields -
                            determines simulation run type
    coords              <Vector{PolyVec}> list of polygon coords to make into floes
    domain              <Domain> model domain 
    hmean               <Float> average floe height
    Δh                  <Float> height range - floes will range in height from
                            hmean ± Δh
    floe_settings       <FloeSettings> settings needed to initialize floes
    rng                 <RNG> random number generator to generate random floe
                            attributes - default uses Xoshiro256++ algorithm
Output:
    floe_arr <StructArray{Floe}> list of floes created from given polygon
    coordinates
"""
function initialize_floe_field(
    ::Type{FT},
    coords,
    domain,
    hmean,
    Δh;
    floe_settings = FloeSettings(min_floe_area = 0.0),
    rng = Xoshiro(),
) where {FT <: AbstractFloat}
    floe_arr = StructArray{Floe{FT}}(undef, 0)
    floe_polys = [LG.Polygon(valid_polyvec!(c)) for c in coords]
    # Remove overlaps with topography
    if !isempty(domain.topography)
        floe_polys = GO.difference(LG.MultiPolygon(floe_polys), LG.MultiPolygon(domain.topography.coords); target = GI.PolygonTrait(),fix_multipoly = nothing)
    end
    # Turn polygons into floes
    for p in floe_polys
        poly_to_floes!(
            FT,
            floe_arr,
            p,
            hmean,
            Δh,
            domain.east.val - domain.west.val;
            floe_settings = floe_settings,
            rng = rng,
        )
    end
    # Warn about floes with area less than minimum floe size
    min_floe_area = floe_settings.min_floe_area > 0 ?
        floe_settings.min_floe_area :
        FT(
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
    )

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
    domain_coords::Vector{<:PolyVec{<:FT}},
    rng,
    min_to_warn::Int;
    max_tries::Int = 10,
) where {FT <: AbstractFloat}
    xpoints = Vector{FT}()
    ypoints = Vector{FT}()
    domain_poly = GI.MultiPolygon(GO.tuples(domain_coords))
    area_frac = GO.area(domain_poly) / reduce(*, scale_fac)
    # Increase the number of points based on availible percent of bounding box
    npoints = ceil(Int, desired_points / area_frac)
    current_points = 0
    tries = 0
    while current_points < desired_points && tries <= max_tries
        x = rand(rng, FT, npoints)
        y = rand(rng, FT, npoints)
        # Check which of the scaled and translated points are within the domain coords
        in_idx = [GO.coveredby(
            (scale_fac[1] * x[i] .+ trans_vec[1], scale_fac[2] * y[i] .+ trans_vec[2]),
            domain_poly
        ) for i in eachindex(x)]
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
    if current_points > 1
        tess_cells = voronoicells(
            xpoints,
            ypoints,
            Rectangle(Point2(0.0, 0.0), Point2(1.0, 1.0)),
            rng = rng
        ).Cells
        # Scale and translate voronoi coordinates
        tcoords = Vector{PolyVec{FT}}(undef, length(tess_cells))
        for i in eachindex(tess_cells)
            perturb_vec = [
                (-1)^rand(rng, 0:1) * rand(rng, FT)*1e-10,
                (-1)^rand(rng, 0:1) * rand(rng, FT)*1e-10,
            ]
            tcoords[i] = [valid_ringvec!([
                Vector(c) .* scale_fac .+
                trans_vec .+ perturb_vec
                for c in tess_cells[i]
            ])]
        end
        return tcoords
    else
        return Vector{PolyVec{FT}}()
    end
end

"""
    initialize_floe_field(
        ::Type{FT},
        nfloes,
        concentrations,
        domain,
        hmean,
        Δh;
        floe_settings,
        rng,
    )

Create a field of floes using Voronoi Tesselation.
Inputs:
    Type{FT}        <AbstractFloat> Type for grid's numberical fields -
                        determines simulation run type
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
    floe_settings       <FloeSettings> settings needed to initialize floes
    rng                 <RNG> random number generator to generate random floe
                            attributes - default uses Xoshiro256++
Output:
    floe_arr <StructArray> list of floes created using Voronoi Tesselation
        of the domain with given concentrations.
"""
function initialize_floe_field(
    ::Type{FT},
    nfloes::Int,
    concentrations,
    domain,
    hmean,
    Δh;
    floe_settings = FloeSettings(min_floe_area = 0.0),
    rng = Xoshiro(),
) where {FT <: AbstractFloat}
    floe_arr = StructArray{Floe{FT}}(undef, 0)
    nfloes_added = 0
    # Split domain into cells with given concentrations
    nrows, ncols = size(concentrations[:, :])
    Lx = domain.east.val - domain.west.val
    Ly = domain.north.val - domain.south.val
    rowlen = Ly / nrows
    collen = Lx / ncols
    # Availible space in whole domain
    topography_list = [LG.Polygon(t) for t in domain.topography.coords]
    open_water_area = (Lx * Ly) - sum(GO.area, topography_list; init = 0.0)
    # Loop over cells
    for j in range(1, ncols)
        for i in range(1, nrows)
            c = concentrations[i, j]
            if c > 0
                c = c > 1 ? 1 : c
                # Grid cell bounds
                xmin = domain.west.val + collen * (j - 1)
                ymin = domain.south.val + rowlen * (i - 1)
                trans_vec = [xmin, ymin]
                # Open water in cell
                open_cell_init = LG.Polygon(rect_coords(xmin, xmin + collen, ymin, ymin + rowlen))
                open_cell = if !isempty(domain.topography)
                    diff_polys(open_cell_init, GI.MultiPolygon(domain.topography.coords); fix_multipoly = nothing)
                else
                    [open_cell_init]
                end
                open_cell_mpoly = GI.MultiPolygon(open_cell)
                open_coords = [GI.coordinates(c) for c in open_cell]
                open_area = sum(GO.area, open_cell; init = 0.0)
                # Generate coords with voronoi tesselation and make into floes
                ncells = ceil(Int, nfloes * open_area / open_water_area / c)
                floe_coords = generate_voronoi_coords(
                    ncells,
                    [collen, rowlen],
                    trans_vec,
                    open_coords,
                    rng,
                    ncells,
                )
                nfloes = length(floe_coords)
                if !isempty(floe_coords)
                    floe_poly_list = [LG.Polygon(c) for c in floe_coords]
                    nfloes = length(floe_poly_list)
                    floe_idx = shuffle(rng, range(1, nfloes))
                    floes_area = FT(0.0)
                    while !isempty(floe_idx) && floes_area/open_area <= c
                        idx = pop!(floe_idx)
                        poly_pieces_list = intersect_polys(floe_poly_list[idx], open_cell_mpoly)
                        for piece in poly_pieces_list
                            n_new_floes = poly_to_floes!(
                                FT,
                                floe_arr,
                                piece,
                                hmean,
                                Δh,
                                domain.east.val - domain.west.val;
                                floe_settings = floe_settings,
                                rng = rng,
                            )
                            floes_area += sum(Iterators.drop(floe_arr.area, nfloes_added))
                            nfloes_added += n_new_floes
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




