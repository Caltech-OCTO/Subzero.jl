"""
Types to hold parameters for simulation's physical processes
"""

"""
    CouplingSettings

Settings needed for coupling within the model.
If coupling_on is true, the model will be coupled with the simulation's ocean
and atmosphere. The Δt determines how many simulation timesteps between
calculating ocean and atmospheric forces on the floes. If calc_ocnτ is true then
the simulation calculates the stress the ice and atmosphere put on the ocean. 
"""
@kwdef struct CouplingSettings
    coupling_on::Bool = true
    Δt::Int = 10
    calc_ocnτ::Bool = false

    function CouplingSettings(coupling_on, Δt, calc_ocnτ)
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
        floe_floe_max_overlap::FT, 
        floe_domain_max_overlap::FT,
    ) where {FT<:AbstractFloat}
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
        new{FT}(
            collisions_on,
            floe_floe_max_overlap,
            floe_domain_max_overlap,
        )
    end
    CollisionSettings(
        collisions_on,
        floe_floe_max_overlap::FT,
        floe_domain_max_overlap::FT,
    ) where {FT<:AbstractFloat} = 
        CollisionSettings{FT}(
            collisions_on,
            floe_floe_max_overlap,
            floe_domain_max_overlap,
        )
end

"""
    FractureSettings{CT<:AbstractFractureCriteria}

Settings needed for fractures within the model. 
If fractures_on is true, fractures will occur, else they will not.
The criteria defines which fracture criteria are used to determine which floes
to fracture. The Δt determines how many simulation timesteps between fracturing
floes. If deform_on is true, then the floe will be deformed around floe
primarily causing the fracture, identified by the largest overlap area on the
most recent set of collisions. Npieces denotes how many pieces to try to split a
fractured floe into - 3 is suggested value.
"""
@kwdef struct FractureSettings{CT<:AbstractFractureCriteria}
    fractures_on::Bool = false
    criteria::CT = NoFracture()
    Δt::Int = 0
    deform_on::Bool = false
    npieces::Int = 1

    function FractureSettings{CT}(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
    ) where {CT <: AbstractFractureCriteria}
        if fractures_on
            if Δt < 0
                @warn "Fracturing can't occur on a multiple of negative timesteps. Turning fracturing off."
                fractures_on = false
            elseif criteria isa NoFracture
                @warn "Fracturing can't occur on with NoFracture criteria. Turning fracturing off."
                fractures_on = false
            elseif npieces < 2
                @warn "Fracturing can't occur on with npieces < 2 as this won't split floe. Turning fracturing off."
                fractures_on = false
            end
        end
        if !fractures_on && deform_on
            @warn "Deformation can't occur on without fracturing. Turning deformation off."
            deform_on = false
        end
        new{CT}(fractures_on, criteria, Δt, deform_on, npieces)
    end
    FractureSettings(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
    ) where {CT <: AbstractFractureCriteria} = 
        FractureSettings{CT}(
            fractures_on,
            criteria,
            Δt,
            deform_on,
            npieces,
        )
end


"""
    struct SimplificationSettings{FT<:AbstractFloat}

Settings needed for simplification within the model.
If dissolve_on is true, floes below the min_floe_area will be removed from the
simulation every timestep. If smooth_vertices_on is true then floe's with
more vertices than max_vertices will be simplified every Δt_smooth timesteps.
"""
@kwdef struct SimplificationSettings{FT<:AbstractFloat}
    dissolve_on::Bool = true
    min_floe_area::FT = 1e4
    smooth_vertices_on::Bool = true
    max_vertices::Int = 30
    Δt_smooth::Int = 20

    function SimplificationSettings{FT}(
        dissolve_on,
        min_floe_area::FT,
        smooth_vertices_on,
        max_vertices,
        Δt_smooth
    ) where {FT<:AbstractFloat}
        if dissolve_on && min_floe_area <= 0
            @warn "If the minimum floe area is less than or equal to 0, no floes will be dissolved. Turning dissolve floes off."
            dissolve_on = false
        end
        if smooth_vertices_on && Δt_smooth < 0
            @warn "Floe smoothing can't occur on a multiple of negative timesteps. Turning floe simplification off."
            smooth_vertices_on = false
        end
        new{FT}(
            dissolve_on,
            min_floe_area,
            smooth_vertices_on,
            max_vertices,
            Δt_smooth,
        )
    end

    SimplificationSettings(
        dissolve_on,
        min_floe_area::FT,
        smooth_vertices_on,
        max_vertices,
        Δt_smooth
    ) where {FT<:AbstractFloat} = 
        SimplificationSettings{FT}(
            dissolve_on,
            min_floe_area,
            smooth_vertices_on,
            max_vertices,
            Δt_smooth
        )
end

"""
    SimplificationSettings(
        grid::AbstractGrid,
        dissolve_on::Bool,
        smooth_vertices_on::Bool,
        max_vertices::Int,
        Δt_smooth::Int,
    )
Determine simplification settings with minimum floe area depending on the size
of the model's grid. 
Inputs:
    grid                <AbstractGrid> model's grid
    dissolve_on         <Bool> true if floes below minimum floe area should be
                            removed from simulation
    smooth_vertices_on  <Bool> true if floes with more vertices than
                            max_vertices should be simplified
    max_vertices        <Int> maximum number of verticies before simplifying
                            floe shape
    Δt_smooth           <Int> number of timesteps between each round of floe
                            vertex simplification
Outputs:
    SimplificationSettings object where the minimum floe area is set to 4 times
    the grid area, divided by 1e4 and the rest of the parameters are user input
"""
function SimplificationSettings(
    grid::AbstractGrid,
    dissolve_on::Bool,
    smooth_vertices_on::Bool,
    max_vertices::Int,
    Δt_smooth::Int,
)
    min_floe_area = 4 * (grid.xg[end] - grid.xg[1]) *
        (grid.yg[end] - grid.yg[1]) / 1e4

    return SimplificationSettings(
        dissolve_on,
        eltype(grid.xg)(min_floe_area),
        smooth_vertices_on,
        max_vertices,
        Δt_smooth,
    ) 
end

# CORNERS::Bool = false           # If true, corners of floes can break
# PACKING::Bool = false           # If true, floe packing is enabled
# RAFTING::Bool = false           # If true, floe rafting is enabled
# RIDGING::Bool = false           # If true, floe ridging is enabled
# WELDING::Bool = false           # If true, floe welding is enabled
#Δtpack::Int = 500               # Timesteps between thermodynamic floe 
                                    # creation