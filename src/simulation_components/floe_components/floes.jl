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
Singular sea ice floe with fields describing current state.
"""
@kwdef mutable struct Floe{FT<:AbstractFloat}
    # Physical Properties -------------------------------------------------
    poly::Polys{FT}         # polygon that represents the floe's shape
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
    stress_accum::Matrix{FT} = zeros(2, 2)
    stress_instant::Matrix{FT} = zeros(2, 2)
    strain::Matrix{FT} = zeros(2, 2)
    damage::FT = 0.0        # damage to the floe (to be used in new stress calculator)
    # Previous values for timestepping  -------------------------------------
    p_dxdt::FT = 0.0        # previous timestep x-velocity
    p_dydt::FT = 0.0        # previous timestep y-velocity
    p_dudt::FT = 0.0        # previous timestep x-acceleration
    p_dvdt::FT = 0.0        # previous timestep y-acceleration
    p_dξdt::FT = 0.0        # previous timestep time derivative of ξ
    p_dαdt::FT = 0.0        # previous timestep angular-velocity
end

const FloeType{FT} = Union{LazyRow{Floe{FT}}, Floe{FT}} where FT

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

Constructor for floe with a polygon
Inputs:
    poly                <Polygon> 
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
    floe = GO.tuples(poly, FT)
    rmholes!(floe)
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
    moment = _calc_moment_inertia(FT, floe, centroid, height; ρi = floe_settings.ρi)
    angles = GO.angles(floe, FT)
    rmax = calc_max_radius(floe, centroid, FT)
    status = Status()
    # Generate Monte Carlo Points
    x_subfloe_points, y_subfloe_points, status = generate_subfloe_points(
        floe_settings.subfloe_point_generator,
        floe,
        centroid,
        area_tot,
        status,
        rng,
    )
    # Generate Stress History
    stress_instant = zeros(FT, 2, 2)

    return Floe{FT}(;
        poly = GO.tuples(floe, FT),
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
        stress_instant = stress_instant,
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
function Floe{FT}(
    coords::PolyVec,
    hmean,
    Δh;
    floe_settings = FloeSettings(),
    rng = Xoshiro(),
    kwargs...,
) where {FT <: AbstractFloat}
    valid_polyvec!(coords)
    rmholes!(coords)
    return Floe{FT}(
        make_polygon(coords),
        hmean,
        Δh;
        floe_settings = floe_settings,
        rng = rng,
        kwargs...,
    ) 

end

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
