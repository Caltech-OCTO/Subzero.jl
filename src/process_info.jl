"""
Types to hold parameters for simulation's physical processes
"""

"""
    CouplingInfo

Parameters needed for coupling within the model.
If coupling_on is true, the model will be coupled with the simulation's ocean and atmosphere.
The Δt determines how many simulation timesteps between calculating ocean
and atmospheric forces on the floes. If calc_ocnτ is true then the simulation
calculates the stress the ice and atmosphere put on the ocean. 
"""
struct CouplingInfo
    coupling_on::Bool
    Δt::Int
    calc_ocnτ::Bool

    function CouplingInfo(coupling_on, Δt, calc_ocnτ)
        if coupling_on && Δt < 0
            @warn "Coupling can't occur on a multiple of negative timesteps. Turning coupling off."
            coupling_on = false
        end
        if !coupling_on && calc_ocnτ
            @warn "Can't calculate stresses on ocean from ice and atmosphere without coupling. Turning calc_ocnτ off."
            calc_ocnτ = false
        end
        new(coupling_on, Δt, calc_ocnτ)
    end
end

"""
    CouplingInfo()

Default coupling information. Forces on floes from ocean and atmosphere are
calculated every 10 timesteps and the stresses on the ocean from the ice and
atmosphere are not calculated. 
"""
CouplingInfo() = CouplingInfo(true, 10, false)

"""
    CollisionInfo{FT<:AbstractFloat}

Parameters needed for collisions within the model. 
If collisions_on is true, collisions will occur, else they will not.
The floe_floe_max_overlap defines the percentage of overlap allowed between
floes before marking them for ridging/rafting.
The floe_domain_max_overlap defines the percentage of overlap allowed between
floes and the domain (collision boundaries and topography) before removing the
floe from the simulation. 
Both floe_floe_max_overlap and floe_domain_max_overlap should be between 0-1 and
if a value below 0 is given or a value greater than 1 is given when collisions_on is
true they will be set to 0 and 1 respectively.
"""
struct CollisionInfo{FT<:AbstractFloat}
    collisions_on::Bool
    floe_floe_max_overlap::FT
    floe_domain_max_overlap::FT

    function CollisionInfo{FT}(collisions_on, floe_floe_max_overlap::FT, 
                                           floe_domain_max_overlap::FT) where {FT<:AbstractFloat}
        if collisions_on
            if floe_floe_max_overlap > 1
                @warn "The maximum collisin overlap between floes can't be greater than 1. Setting to 1."
                floe_floe_max_overlap = FT(1)
            elseif floe_floe_max_overlap < 0
                @warn "The maximum collisin overlap between floes can't be less than 0. Setting to 0."
                floe_floe_max_overlap = FT(0)
            end

            if floe_domain_max_overlap > 1
                @warn "The maximum collisin overlap between floes and the domain can't be greater than 1. Setting to 1."
                floe_domain_max_overlap = FT(1)
            elseif floe_domain_max_overlap < 0
                @warn "The maximum collisin overlap between floes and the domain can't be less than 0. Setting to 0."
                floe_domain_max_overlap = FT(0)
            end
        end
        new{FT}(collisions_on, floe_floe_max_overlap, floe_domain_max_overlap)
    end
    CollisionInfo(collisions_on, floe_floe_max_overlap::FT, floe_domain_max_overlap::FT) where {FT<:AbstractFloat} = 
        CollisionInfo{FT}(collisions_on, floe_floe_max_overlap, floe_domain_max_overlap)
end

"""
    CollisionInfo()

Default collision information. Collisions are turned on and the maximum
floe-floe overlap is 55% (0.55) and the maximum floe-domain overlap is 75% (0.75).
"""
CollisionInfo() = CollisionInfo(true, 0.55, 0.75)

"""
    FractureInfo{CT<:AbstractFractureCriteria}

Parameters needed for fractures within the model. 
If fractures_on is true, fractures will occur, else they will not.
The criteria defines which fracture criteria are used to determine which floes
to fracture. The Δt determines how many simulation timesteps between fracturing
floes. If deform_on is true, then the floe will be deformed around floe
primarily causing the fracture, identified by the largest overlap area on the
most recent set of collisions. 
"""
struct FractureInfo{CT<:AbstractFractureCriteria}
    fractures_on::Bool
    criteria::CT
    Δt::Int
    deform_on::Bool

    function FractureInfo{CT}(fractures_on, criteria::CT, Δt, deform_on) where {CT <: AbstractFractureCriteria}
        if fractures_on
            if Δt < 0
                @warn "Fracturing can't occur on a multiple of negative timesteps. Turning fracturing off."
                fractures_on = false
            elseif criteria isa NoFracture
                @warn "Fracturing can't occur on with NoFracture criteria. Turning fracturing off."
                fractures_on = false
            end
        end
        if !fractures_on && deform_on
            @warn "Deformation can't occur on without fracturing. Turning deformation off."
            deform_on = false
        end
        new{CT}(fractures_on, criteria, Δt, deform_on)
    end
    FractureInfo(fractures_on, criteria::CT, Δt, deform_on) where {CT <: AbstractFractureCriteria} = 
        FractureInfo{CT}(fractures_on, criteria, Δt, deform_on)

end

"""
    FractureInfo()

Default fracture information. Fractures are turned off with fractures_on set to
false. The NoFracture criteria is provided, Δt is set to 0, and deform_on is
also false, although none of these will be used. 
"""
FractureInfo() = FractureInfo(false, NoFracture(), 0, false)

# CORNERS::Bool = false           # If true, corners of floes can break
# KEEPMIN::Bool = false           # If true, retain small floes that would 
#                                 # normally "dissolve"
# PACKING::Bool = false           # If true, floe packing is enabled
# RAFTING::Bool = false           # If true, floe rafting is enabled
# RIDGING::Bool = false           # If true, floe ridging is enabled
# WELDING::Bool = false           # If true, floe welding is enabled
#Δtsimp::Int = 20                # Timesteps between floe simplification
#Δtpack::Int = 500               # Timesteps between thermodynamic floe 
                                    # creation