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
height. If smooth_vertices_on is true then floe's with more vertices than
max_vertices will be simplified every Δt_smooth timesteps.

A reasonable formula for minimum floe area is the following:
min_floe_area = 4(grid.xg[end] - grid.xg[1]) * (grid.yg[end] - grid.yg[1]) / 1e4
"""
@kwdef struct SimplificationSettings{FT<:AbstractFloat}
    min_floe_area::FT = 1e6
    min_floe_height::FT = 0.1
    max_floe_height::FT = 10.0
    smooth_vertices_on::Bool = true
    max_vertices::Int = 30
    tol::FT = 100.0  # Douglas–Peucker algorithm tolerance in (m)
    Δt_smooth::Int = 20

    function SimplificationSettings{FT}(
        min_floe_area,
        min_floe_height,
        max_floe_height,
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
        new{FT}(
            min_floe_area,
            min_floe_height,
            max_floe_height,
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

Settings needed for ridging and rafting within the model.
"""
@kwdef struct RidgeRaftSettings{FT<:AbstractFloat}
    ridge_raft_on::Bool = false
    Δt::Int = 0
    ridge_probability::FT = 0.95
    raft_probability::FT = 0.95
    min_overlap::FT = 500
    min_ridge_height::FT = 0.2
    max_ridge_height::FT = 5.0
    max_raft_height::FT = 0.25
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

# CORNERS::Bool = false           # If true, corners of floes can break
# PACKING::Bool = false           # If true, floe packing is enabled
# WELDING::Bool = false           # If true, floe welding is enabled
#Δtpack::Int = 500               # Timesteps between thermodynamic floe 
                                    # creation