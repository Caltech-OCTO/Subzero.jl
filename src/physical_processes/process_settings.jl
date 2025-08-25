export CouplingSettings, CollisionSettings, FractureSettings
export SimplificationSettings, RidgeRaftSettings, WeldSettings

"""
    CouplingSettings

Settings needed for coupling within the model.

## _Fields_
- `coupling_on::Bool`: if true, the model will be coupled with the simulation's ocean and atmosphere (Default = true)
- `Δt::Int`: determines how many simulation timesteps between calculating ocean and atmospheric forces on the floes (Default = 10)
- `Δd::Int`: Δd number of buffer grid cells on each side of floe for interpolation using sub-floe points - see [`AbstractSubFloePointsGenerator`](@ref) (Default = 1)
- `two_way_coupling_on::Bool`: if true, then the simulation calculates the stress the ice/atmosphere put on the ocean so that the user could couple to Oceananigans (Default = False) 

!!! note
    `two_way_coupling_on` does NOT set up the coupling with Oceananigans. The user still must do that. This simply calcualtes the fields that are needed for coupling so that they can be passed.

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `CouplingSettings` 
```jldoctest
julia> CouplingSettings()
CouplingSettings(true, 10, 1, false)
```

- Creating `CouplingSettings` with a  `two_way_coupling_on` and a larger interpolation buffer.
```jldoctest
julia> CouplingSettings(Δd = 2, two_way_coupling_on = true)
CouplingSettings(true, 10, 2, true)
```
"""
@kwdef struct CouplingSettings
    coupling_on::Bool = true
    Δt::Int = 10
    Δd::Int = 1
    two_way_coupling_on::Bool = false

    function CouplingSettings(
        coupling_on,
        Δt,
        Δd,
        two_way_coupling_on,
    )
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
        new(
            coupling_on,
            Δt,
            Δd,
            two_way_coupling_on,
        )
    end
end

# See below for docs
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
    CollisionSettings(
        collisions_on,
        floe_floe_max_overlap,
        floe_domain_max_overlap,
    ) = CollisionSettings{Float64}(
        collisions_on,
        floe_floe_max_overlap,
        floe_domain_max_overlap,
    )
end

"""
    CollisionSettings{FT}

Settings needed for collisions within the model. 

## _Fields_
- `collisions_on::Bool`: if true, collisions between floes and the domain elements will occur (Default = true)
- `floe_floe_max_overlap::FT`: defines the percentage of overlap allowed between floes before marking them for ridging/rafting (Default = 0.55)
- `floe_domain_max_overlap::FT`: defines the percentage of overlap allowed between floes and [`DomainElements`](@ref) before removing the floes from the simulation. 

Both `floe_floe_max_overlap` and `floe_domain_max_overlap` should be between 0-1 and if a value < 0 is given or a value > 1 is given when `collisions_on` is true they
will be set to 0 and 1 respectively.

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `CollisionSettings` 
```jldoctest
julia> CollisionSettings()
CollisionSettings{Float64}(true, 0.55, 0.75)
```

- Creating a Float32 `CollisionSettings` with a higher `floe_floe_max_overlap`
```jldoctest
julia> CollisionSettings(Float32;  floe_floe_max_overlap = 0.65)
CollisionSettings{Float32}(true, 0.65f0, 0.75f0)
```
"""
CollisionSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    CollisionSettings{FT}(;kwargs...)

"""
    FractureSettings{CT}

Settings needed for fractures within the model. 

## _Fields_
- `fractures_on::Bool`: if true, floes will fracture, else they will not (Default = false)
- `criteria::CT`: defines which fracture criteria (subtype of [`AbstractFractureCriteria`](@ref)) are used to determine which floes to fracture
- `Δt::Int`: determines how many simulation timesteps between attempting floe fracture (Default = 0)
- `deform_on::Bool`: if true, then the floe will be deformed around floe primarily causing the fracture, identified by the largest overlap area on the most recent set of collisions (Default = False)
- `npieces::Int`: how many pieces to try to split a fractured floe into (Default = 3)

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `FractureSettings` 
```jldoctest
julia> FractureSettings()
FractureSettings{NoFracture}
  ⊢ fractures_on = false
  ⊢ criteria = NoFracture
  ⊢ Δt = 0
  ⊢ deform_on = false
  ⊢ npieces = 3
```

- Creating a Float32 `FractureSettings` with fractures on
```jldoctest
julia> mohrs_cone_criteria = MohrsCone(Float32);

julia> FractureSettings(; fractures_on = true, criteria = mohrs_cone_criteria)
FractureSettings{MohrsCone{Float32}}
  ⊢ fractures_on = true
  ⊢ criteria = MohrsCone{Float32}
  ⊢ Δt = 0
  ⊢ deform_on = false
  ⊢ npieces = 3
```
"""
@kwdef struct FractureSettings{CT<:AbstractFractureCriteria}
    fractures_on::Bool = false
    criteria::CT = NoFracture()
    Δt::Int = 0
    deform_on::Bool = false
    npieces::Int = 3

    function FractureSettings{CT}(
        fractures_on,
        criteria::CT,
        Δt,
        deform_on,
        npieces,
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

# Pretty printing for FractureSettings showing key dimensions
function Base.show(io::IO, fracture_settings::FractureSettings{CT}) where {CT}
    overall_summary = "FractureSettings{$CT}"
    consts_summary = ""
    for name in fieldnames(FractureSettings)
        if name == :criteria
            consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(typeof(getfield(fracture_settings, name))) * "\n")
        else
            consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(getfield(fracture_settings, name)) * "\n")
        end
    end
    print(io, overall_summary, "\n", consts_summary)
end

# See below of documentation
@kwdef struct SimplificationSettings{FT<:AbstractFloat}
    smooth_vertices_on::Bool = true
    max_vertices::Int = 30
    tol::FT = 100.0
    Δt_smooth::Int = 20

    function SimplificationSettings{FT}(
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
            smooth_vertices_on,
            max_vertices,
            tol,
            Δt_smooth,
        )
    end

    SimplificationSettings(
        smooth_vertices_on,
        max_vertices,
        tol,
        Δt_smooth
    ) = SimplificationSettings{Float64}(
        smooth_vertices_on,
        max_vertices,
        tol,
        Δt_smooth
    )
end

"""
    SimplificationSettings{FT}

Settings needed for floe simplification within the simulation.

## _Fields_

- `smooth_vertices_on::Bool`: if true then floe's with more vertices than `max_vertices` will be simplified every `Δt_smooth` timesteps (Default = True)
- `max_vertices::Int`: total number of vertices a floe can have without triggering smoothing if `smooth_vertices_on = True` (Default = 30)
- `tol::FT`: tolerance in meters of the the Douglas–Peucker algorithm, which is used to smooth edges of polygons (Default = 100)
- `Δt_smooth::Int`: determines how many simulation timesteps between attempting floe smoothing (Default = 20)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `SimplificationSettings` 
```jldoctest
julia> SimplificationSettings()
SimplificationSettings{Float64}(true, 30, 100.0, 20)
```

- Creating a Float32 `SimplificationSettings` with `smooth_vertices_on` off
```jldoctest
julia> SimplificationSettings(Float32; smooth_vertices_on = false)
SimplificationSettings{Float32}(false, 30, 100.0f0, 20)
```
"""
SimplificationSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    SimplificationSettings{FT}(;kwargs...)

# See documentation below
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

    RidgeRaftSettings(
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
        domain_gain_probability
    ) = RidgeRaftSettings{Float64}(
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

"""
    RidgeRaftSettings{FT}

Settings needed for ridging and rafting within the simulation.

## _Fields_

- `ridge_raft_on::Bool`: if true, ridging and rafting should be turned on in the simulation (Default = False)
- `Δt`: determines how many simulation timesteps between attempting ridging and rafting (Default = 0)
- `ridge_probability::FT`: the probability a floe ridges with another floe/domain if it meets all other criteria
- `raft_probability::FT`: the probability a floe rafts with another floe/domain if it meets all other criteria
- `min_overlap_frac::FT`: the minimum overlap area fraction between a floe and another floe/domain for that floe to ridge or raft
- `min_ridge_height::FT`: the minimum floe height to ridge with a floe/domain
- `max_floe_ridge_height::FT`: the maximum floe height to ridge with another floe
- `max_domain_rdige_height::FT`: maximum floe height to ridge with a domain element
- `max_floe_raft_height::FT`: maximum floe height to raft with another floe
- `max_domain_raft_height::FT`: maximum floe height to raft with a domain element
- `domain_gain_probability::FT`: the probalility that a floe that rafts with a domain element keeps all of its mass (0) or 
  if that mass is removed and lost to the domain element (1).

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `RidgeRaftSettings` 
```jldoctest
julia> RidgeRaftSettings()
RidgeRaftSettings{Float64}
  ⊢ ridge_raft_on = false
  ⊢ Δt = 0
  ⊢ ridge_probability = 0.95
  ⊢ raft_probability = 0.95
  ⊢ min_overlap_frac = 0.01
  ⊢ min_ridge_height = 0.2
  ⊢ max_floe_ridge_height = 5.0
  ⊢ max_domain_ridge_height = 1.25
  ⊢ max_floe_raft_height = 0.25
  ⊢ max_domain_raft_height = 0.25
  ⊢ domain_gain_probability = 1.0
```

- Creating a Float32 `RidgeRaftSettings` with ridging and rafting turned on  
```jldoctest
julia> RidgeRaftSettings(Float32; ridge_raft_on = true)
RidgeRaftSettings{Float32}
  ⊢ ridge_raft_on = true
  ⊢ Δt = 0
  ⊢ ridge_probability = 0.95
  ⊢ raft_probability = 0.95
  ⊢ min_overlap_frac = 0.01
  ⊢ min_ridge_height = 0.2
  ⊢ max_floe_ridge_height = 5.0
  ⊢ max_domain_ridge_height = 1.25
  ⊢ max_floe_raft_height = 0.25
  ⊢ max_domain_raft_height = 0.25
  ⊢ domain_gain_probability = 1.0
```
"""
RidgeRaftSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    RidgeRaftSettings{FT}(;kwargs...)

# Pretty printing for RidgeRaftSettings showing key dimensions
function Base.show(io::IO, ridge_raft_settings::RidgeRaftSettings{FT}) where {FT}
    overall_summary = "RidgeRaftSettings{$FT}"
    consts_summary = ""
    for name in fieldnames(RidgeRaftSettings)
        consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(getfield(ridge_raft_settings, name)) * "\n")
    end
    print(io, overall_summary, "\n", consts_summary)
end

# See documentation below
@kwdef struct WeldSettings{FT<:AbstractFloat}
    weld_on::Bool = false
    Δts::Vector{Int} = Vector{Int}()
    Nxs::Vector{Int} = Vector{Int}()
    Nys::Vector{Int} = Vector{Int}()
    min_weld_area::FT = 1e6
    max_weld_area::FT = 2e9
    welding_coeff::FT = 150
    function WeldSettings{FT}(
        weld_on,
        Δts,
        Nxs,
        Nys,
        min_weld_area,
        max_weld_area,
        welding_coeff,
    ) where {FT<:AbstractFloat}
        if weld_on && (isempty(Δts) || any(Δts .≤ 0))
            @warn "Welding can't occur without any given timesteps or with \
            negative timesteps. Turning welding off."
            weld_on = false
        elseif any(Nxs .< 1) || any(Nys .< 1)
            @warn "Can't split the grid into less than one row or column. \
            Turning welding off." 
            weld_on = false
        elseif !(length(Δts) == length(Nxs) == length(Nys))
            @warn "Length of timestep multiple list (Δts) must match length of \
            grid split lists Nxs and Nys. Turning welding off." 
            weld_on = false
        end
        # Sort by largest to smallest timestep multiples
        order = reverse!(sortperm(Δts))
        Δts .= Δts[order]
        Nxs .= Nxs[order]
        Nys .= Nys[order]
        new{FT}(
            weld_on,
            Δts,
            Nxs,
            Nys,
            min_weld_area,
            max_weld_area,
            welding_coeff,
        )
    end

    WeldSettings(
        weld_on,
        Δts,
        Nxs,
        Nys,
        min_weld_area,
        max_weld_area,
        welding_coeff,
    ) = WeldSettings{Float64}(
        weld_on,
        Δts,
        Nxs,
        Nys,
        min_weld_area,
        max_weld_area,
        welding_coeff,
    )
end

"""
    WeldSettings{FT}

Settings needed for welding within the simulation.

## _Fields_

- `weld_on::Bool`: a boolean flag for if welding should be turned on in the simulation (Default = False)
- `Δts::Vector{Int}`: a list of multiples of timesteps during which welding code will run,
    welding will be run at multiples of all elements, each with domain split (Default = Int[])
    into corresponding Nx and Ny values
- `Nxs::Vector{Int}`: a list of number of x-directional bins to split the domain into at
    corresponding timesteps (Default = Int[])
- `Nys::Vector{Int}`: a list of number of x-directional bins to split the domain into at
    corresponding timesteps (Default = Int[])
- `min_weld_area::FT`: minimum area a weld can create for two floes to weld (Default = 1e6)
- `max_weld_area::FT`: maximum area a weld can create for two floes to weld (Default = 2e9)
- `welding_coeff::FT`: non-dimensional parameter, multiplied by ratio of overlap
    between two floes to original floe area to determin probability that a floe
    will merge. The larger this is, the more likely floes are to weld.
    Probability with 5% overlap is `welding_coeff * (0.05) > rand()` (Default = 150)

## _Positional arguments_
- $FT_DEF

## _Keyword arguments_
- Each of the above fields is an optional keyword argument.

## _Examples_
- Creating default `WeldSettings` 
```jldoctest
julia> WeldSettings()
WeldSettings{Float64}
  ⊢ weld_on = false
  ⊢ Δts = Int64[]
  ⊢ Nxs = Int64[]
  ⊢ Nys = Int64[]
  ⊢ min_weld_area = 1.0e6
  ⊢ max_weld_area = 2.0e9
  ⊢ welding_coeff = 150.0
```
- Creating a Float32 `WeldSettings` with welding turned on 
```jldoctest
julia> WeldSettings(Float32; weld_on = true, Δts = [100, 200], Nxs = [2, 1], Nys = [1, 2])
WeldSettings{Float32}
  ⊢ weld_on = true
  ⊢ Δts = [200, 100]
  ⊢ Nxs = [1, 2]
  ⊢ Nys = [2, 1]
  ⊢ min_weld_area = 1.0e6
  ⊢ max_weld_area = 2.0e9
  ⊢ welding_coeff = 150.0
```
"""
WeldSettings(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    WeldSettings{FT}(; kwargs...)

# Pretty printing for WeldSettings showing key dimensions
function Base.show(io::IO, weld_settings::WeldSettings{FT}) where {FT}
    overall_summary = "WeldSettings{$FT}"
    consts_summary = ""
    for name in fieldnames(WeldSettings)
        consts_summary *= "  ⊢" * " " * (string(name) * " = " * string(getfield(weld_settings, name)) * "\n")
    end
    print(io, overall_summary, "\n", consts_summary)
end


# CORNERS::Bool = false           # If true, corners of floes can break
# PACKING::Bool = false           # If true, floe packing is enabled