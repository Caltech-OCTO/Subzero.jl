"""
Types to hold parameters for simulation's physical process settings
"""

"""
    CouplingSettings

Settings needed for coupling within the model.
If coupling_on is true, the model will be coupled with the simulation's ocean
and atmosphere. The Δt determines how many simulation timesteps between
calculating ocean and atmospheric forces on the floes. Δd number of buffer grid
cells on each side of floe for monte carlo interpolation and mc_n is the number
of monte carlo points to attempt to generage for each floe. If two_way_coupling_on is
true then the simulation calculates the stress the ice/atmosphere put on the
ocean. 
"""
@kwdef struct CouplingSettings{GT <: AbstractSubFloePointsGenerator}
    coupling_on::Bool = true
    Δt::Int = 10
    Δd::Int = 1
    subfloe_point_generator::GT = MonteCarloPointsGenerator()
    two_way_coupling_on::Bool = false

    function CouplingSettings{GT}(
        coupling_on,
        Δt,
        Δd,
        subfloe_point_generator,
        two_way_coupling_on,
    ) where {GT <: AbstractSubFloePointsGenerator}
        if coupling_on && Δt < 0
            @warn "Coupling can't occur on a multiple of negative timesteps. \
                Turning coupling off."
            coupling_on = false
        end
        if !coupling_on && two_way_coupling_on
            @warn "Can't calculate stresses on ocean from ice and atmosphere \
                without coupling. Turning two_way_coupling_on off."
            two_way_coupling_on = false
        end
        if Δd < 0
            @warn "Can't complete interpolation of ocean and atmosphere forces \
                with a buffer of less than 0 grid cells. Setting Δd = 0."
            Δd = 0
        end
        new{GT}(
            coupling_on,
            Δt,
            Δd,
            subfloe_point_generator,
            two_way_coupling_on,
        )
    end

    CouplingSettings(
        coupling_on,
        Δt,
        Δd,
        subfloe_point_generator::GT,
        two_way_coupling_on,
    ) where {GT <: AbstractSubFloePointsGenerator} = 
        CouplingSettings{GT}(
            coupling_on,
            Δt,
            Δd,
            subfloe_point_generator::GT,
            two_way_coupling_on,
        )
end

"""
    CollisionSettings{FT<:AbstractFloat}

Settings needed for collisions within the model. 
If collisions_on is true, collisions will occur, else they will not.
The floe_floe_max_overlap defines the percentage of overlap allowed between
floes before marking them for ridging/rafting.
The floe_domain_max_overlap defines the percentage of overlap allowed between
floes and the domain (collision boundaries and topography) before removing the
floe from the simulation. 
Both floe_floe_max_overlap and floe_domain_max_overlap should be between 0-1 and
if a value < 0 is given or a value > 1 is given when collisions_on is true they
will be set to 0 and 1 respectively.
"""
@kwdef struct CollisionSettings{FT<:AbstractFloat}
    collisions_on::Bool = true
    floe_floe_max_overlap::FT = 0.55
    floe_domain_max_overlap::FT = 0.75

    function CollisionSettings{FT}(
        collisions_on,
        floe_floe_max_overlap, 
        floe_domain_max_overlap,
    ) where {FT<:AbstractFloat}
        if collisions_on
            if floe_floe_max_overlap > 1
                @warn "The maximum collisin overlap between floes can't be \
                    greater than 1. Setting to 1."
                floe_floe_max_overlap = FT(1)
            elseif floe_floe_max_overlap < 0
                @warn "The maximum collisin overlap between floes can't be \
                    less than 0. Setting to 0."
                floe_floe_max_overlap = FT(0)
            end

            if floe_domain_max_overlap > 1
                @warn "The maximum collisin overlap between floes and the \
                    domain can't be greater than 1. Setting to 1."
                floe_domain_max_overlap = FT(1)
            elseif floe_domain_max_overlap < 0
                @warn "The maximum collisin overlap between floes and the \
                    domain can't be less than 0. Setting to 0."
                floe_domain_max_overlap = FT(0)
            end
        end
        new{FT}(
            collisions_on,
            floe_floe_max_overlap,
            floe_domain_max_overlap,
        )
    end
end

"""
    CollisionSettings(::Type{FT}, args...)

A float type FT can be provided as the first argument of any CollisionSettings
constructor. A CollisionSettings of type FT will be created by passing all
other arguments to the correct constructor. 
"""
CollisionSettings(::Type{FT}, args...) where {FT <: AbstractFloat} =
    CollisionSettings{FT}(args...)

"""
    CollisionSettings(args...)

If a type isn't specified, CollisionSettings will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
CollisionSettings(args...) = CollisionSettings{Float64}(args...)

"""
    FractureSettings{CT<:AbstractFractureCriteria}

Settings needed for fractures within the model. 
If fractures_on is true, fractures will occur, else they will not.
The criteria defines which fracture criteria are used to determine which floes
to fracture. The Δt determines how many simulation timesteps between fracturing
floes. If deform_on is true, then the floe will be deformed around floe
primarily causing the fracture, identified by the largest overlap area on the
most recent set of collisions. Npieces denotes how many pieces to try to split a
fractured floe into - 3 is suggested value. nhistory is the number of previous
stress values to hold in each floe's stress history field. The higher the
number, the longer it will take for a floe to fracture.
"""
@kwdef struct FractureSettings{CT<:AbstractFractureCriteria}
    fractures_on::Bool = false
    criteria::CT = NoFracture()
    Δt::Int = 0
    deform_on::Bool = false
    npieces::Int = 3
    nhistory::Int = 1000

    function FractureSettings{CT}(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
        nhistory,
    ) where {CT <: AbstractFractureCriteria}
        if fractures_on
            if Δt < 0
                @warn "Fracturing can't occur with negative timesteps. Turning \
                    fracturing off."
                fractures_on = false
            elseif criteria isa NoFracture
                @warn "Fracturing can't occur on with NoFracture criteria. \
                    Turning fracturing off."
                fractures_on = false
            elseif npieces < 2
                @warn "Fracturing can't occur on with npieces < 2 as this \
                    won't split floe. Turning fracturing off."
                fractures_on = false
            end
        end
        if !fractures_on && deform_on
            @warn "Deformation can't occur on without fracturing. Turning \
                deformation off."
            deform_on = false
        end
        if nhistory < 1
            @warn "Stress history must have at least one element. Setting \
                nhistory = 1."
            nhistory = 1
        end
        new{CT}(fractures_on, criteria, Δt, deform_on, npieces, nhistory)
    end
    FractureSettings(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
        nhistory,
    ) where {CT <: AbstractFractureCriteria} = 
        FractureSettings{CT}(
            fractures_on,
            criteria,
            Δt,
            deform_on,
            npieces,
            nhistory,
        )
end


"""
    struct SimplificationSettings{FT<:AbstractFloat}

Settings needed for simplification within the model.
Floes below the min_floe_area and below min_floe_height will be removed from the
simulation every timestep. Floes above max_floe_height will not be able to gain
height. Floes with aspect ratio less than min_aspect_ratio will be removed.
If smooth_vertices_on is true then floe's with more vertices than
max_vertices will be simplified every Δt_smooth timesteps.

A reasonable formula for minimum floe area is the following:
min_floe_area = 4(grid.xg[end] - grid.xg[1]) * (grid.yg[end] - grid.yg[1]) / 1e4
"""
@kwdef struct SimplificationSettings{FT<:AbstractFloat}
    min_floe_area::FT = 1e6
    min_floe_height::FT = 0.1
    max_floe_height::FT = 10.0
    min_aspect_ratio::FT = 0.05
    smooth_vertices_on::Bool = true
    max_vertices::Int = 30
    tol::FT = 100.0  # Douglas–Peucker algorithm tolerance in (m)
    Δt_smooth::Int = 20

    function SimplificationSettings{FT}(
        min_floe_area,
        min_floe_height,
        max_floe_height,
        min_aspect_ratio,
        smooth_vertices_on,
        max_vertices,
        tol,
        Δt_smooth
    ) where {FT<:AbstractFloat}
        if smooth_vertices_on && Δt_smooth < 0
            @warn "Floe smoothing can't occur on a multiple of negative \
                timesteps. Turning floe simplification off."
            smooth_vertices_on = false
        end
        if min_aspect_ratio > 1
            @warn "Floes can't have an aspect ratio greater than 1. Setting \
            max_aspect_ratio to 1. However, this means only with equal x and y \
            maximum extents can exist. It is suggested this is set much lower \
            (e.g. 0.05)."
            min_aspect_ratio = FT(1)
        elseif min_aspect_ratio < 0
            @warn "Floes can't have an aspect ratio less than 0. Setting \
            max_aspect_ratio to 0. This means that floes of any aspect ratio \
            can exist throughout the simulation."
            min_aspect_ratio = FT(0)
        end
        new{FT}(
            min_floe_area,
            min_floe_height,
            max_floe_height,
            min_aspect_ratio,
            smooth_vertices_on,
            max_vertices,
            tol,
            Δt_smooth,
        )
    end
end

"""
    SimplificationSettings(::Type{FT}, args...)

A float type FT can be provided as the first argument of any
SimplificationSettings constructor. A SimplificationSettings of type FT will be
created by passing all other arguments to the correct constructor. 
"""
SimplificationSettings(::Type{FT}, args...) where {FT <: AbstractFloat} =
    SimplificationSettings{FT}(args...)

"""
    SimplificationSettings(args...)

If a type isn't specified, SimplificationSettings will be of type Float64 and
the correct constructor will be called with all other arguments.
"""
SimplificationSettings(args...) = SimplificationSettings{Float64}(args...)

"""
    RidgeRaftSettings

Settings needed for ridging and rafting within the model. These have the
following meanings:
- ridge_probability: the probability a floe ridges with another floe/domain if
    it meets all other criteria
- raft_probability: the probability a floe rafts with another floe/domain if it
    meets all other criteria
- min_overlap_frac: the minimum overlap area fraction between a floe and another
    floe/domain for that floe to ridge or raft
- min_ridge_height: the minimum floe height to ridge with a floe/domain
- max_floe_ridge_height: the maximum floe height to ridge with another
    floe
- max_domain_rdige_height: maximum floe height to ridge with a domain element
- max_floe_raft_height: maximum floe height to raft with another floe
- max_domain_raft_height: maximum floe height to raft with a domain element
- domain_gain_probability: the probalility that a floe that rafts with a domain
    element keeps all of its mass (0) or if that mass is removed and lost to the
    domain element (1).
"""
@kwdef struct RidgeRaftSettings{FT<:AbstractFloat}
    ridge_raft_on::Bool = false
    Δt::Int = 0
    ridge_probability::FT = 0.95
    raft_probability::FT = 0.95
    min_overlap_frac::FT = 0.01
    min_ridge_height::FT = 0.2
    max_floe_ridge_height::FT = 5.0
    max_domain_ridge_height::FT = 1.25
    max_floe_raft_height::FT = 0.25
    max_domain_raft_height::FT = 0.25
    domain_gain_probability::FT = 1.0

    function RidgeRaftSettings{FT}(
        ridge_raft_on,
        Δt,
        ridge_probability,
        raft_probability,
        min_overlap_frac,
        min_ridge_height,
        max_floe_ridge_height,
        max_domain_ridge_height,
        max_floe_raft_height,
        max_domain_raft_height,
        domain_gain_probability,
    ) where {FT<:AbstractFloat}
        if ridge_raft_on && Δt < 0
            @warn "Ridging and rafting can't occur on a multiple of negative \
                timesteps. Turning ridging and rafting off."
            ridge_raft_on = false
        end
        if ridge_probability > 1
            @warn "Floes can't have a greater ridge probability than 1. \
            Setting ridge probability to 1."
            ridge_probability = FT(1)
        elseif ridge_probability < 0
            @warn "Floes can't have a smaller ridge probability than 1. \
            Setting ridge probability to 0."
            ridge_probability = FT(0)
        end
        if raft_probability > 1
            @warn "Floes can't have a greater raft probability than 1. \
            Setting ridge probability to 1."
            raft_probability = FT(1)
        elseif raft_probability < 0
            @warn "Floes can't have a smaller raft probability than 1. \
            Setting ridge probability to 0."
            raft_probability = FT(0)
        end
        if min_overlap_frac > 1
            @warn "Floes can't overlap more than 100%, so min_overlap_frac \
            can't exceed 1. Setting min_overlap_frac to 1."
            min_overlap_frac = FT(1)
        elseif min_overlap_frac < 0
            @warn "Floes can't overlap less than 0%, so min_overlap_frac \
            can't be less than 0. Setting min_overlap_frac to 0."
            min_overlap_frac = FT(0)
        end
        if domain_gain_probability > 1
            @warn "Floes can't have a greater domain_gain_probability than 1. \
            Setting domain_gain_probability to 1."
            domain_gain_probability = FT(1)
        elseif domain_gain_probability < 0
            @warn "Floes can't have a smaller domain_gain_probability than 0. \
            Setting domain_gain_probability to 0."
            domain_gain_probability = FT(0)
        end
        new{FT}(
            ridge_raft_on,
            Δt,
            ridge_probability,
            raft_probability,
            min_overlap_frac,
            min_ridge_height,
            max_floe_ridge_height,
            max_domain_ridge_height,
            max_floe_raft_height,
            max_domain_raft_height,
            domain_gain_probability,
        )
    end
end

"""
    RidgeRaftSettings(::Type{FT}, args...)

A float type FT can be provided as the first argument of any RidgeRaftSettings
constructor. A RidgeRaftSettings of type FT will be created by passing all other
arguments to the correct constructor. 
"""
RidgeRaftSettings(::Type{FT}, args...) where {FT <: AbstractFloat} =
    RidgeRaftSettings{FT}(args...)

"""
    RidgeRaftSettings(args...)

If a type isn't specified, RidgeRaftSettings will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
RidgeRaftSettings(args...) = RidgeRaftSettings{Float64}(args...)


@kwdef struct WeldSettings{FT<:AbstractFloat}
    max_welded_floe_area::FT
    min_weld_area::FT
    max_weld_area::FT 
    welding_coeff::FT = 150
    Δts::Vector{FT} = [0]
    Nxs::Vector{FT} = [1]
    Nys::Vector{FT} = [1]
end

"""
    WeldSettings(::Type{FT}, args...)

A float type FT can be provided as the first argument of any WeldSettings
constructor. A WeldSettings of type FT will be created by passing all other
arguments to the correct constructor. 
"""
WeldSettings(::Type{FT}, args...) where {FT <: AbstractFloat} =
    WeldSettings{FT}(args...)

"""
    WeldSettings(args...)

If a type isn't specified, WeldSettings will be of type Float64 and the
correct constructor will be called with all other arguments.
"""
WeldSettings(args...) = WeldSettings{Float64}(args...)


# CORNERS::Bool = false           # If true, corners of floes can break
# PACKING::Bool = false           # If true, floe packing is enabled