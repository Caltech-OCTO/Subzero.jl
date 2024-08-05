# # DecayAreaScaledStressCalculator and DamageStressCalculator

export AbstractStressCalculator, DecayAreaScaledStressCalculator, DamageStressCalculator

#= 
## What is a stress calculator? 

A stress calculator defines how the stress for a floe is calculated from its interactions
and how that stress is used to determine fracture. 

The `AbstractStressCalculator` is the abstract type, with `DecayAreaScaledStressCalculator`
currently implemented and the infrastructure for `DamageStressCalculator` setup for future
implementation.

Each of these concrete types represent different methods of calculating a floe's stress. A
`FloeSettings` struct requires a concrete subtype of `AbstractStressCalculator`. The
functions for calculating stress are then dispatched off of this calculator type. The
only function currently is: `_update_stress!`. This function need to be extended for
any new concrete subtype of `AbstractStressCalculator`.
=#

const STRESS_CALCULATOR_DEF = ": calculator to determine floe's stress calculations from floe interactions."

"""
    abstract type AbstractStressCalculator{FT <: AbstractFloat}

Abstract super type for stress calculators, which calculate a floe's stress from floe
interactions. The subtypes serve as dispatch types for the following two calculation
methods.

## API
The following method must be implemented for all subtypes:
- `_update_stress!(stress_calculator::AbstractStressCalculator{FT}, curr_stress::Matrix{FT} , floe::FloeType{FT})`

`_update_stress!` is called in the `calc_stress!` function and takes the stress
calculator, the `floe`'s instantatious stress at the current timestep, and the `floe` itself
and updates the `floe`'s `stress_accum` field, which is used when determining floe fracture
based off of stress. Within the function, other floe fields can be updated as needed.
"""
abstract type AbstractStressCalculator{FT <: AbstractFloat}  end

const DECAY_ARG_WARNING = "λ must be between 0 and 1. Resetting to default value of 0.2."

"""
    DecayAreaScaledStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}

Type of AbstractStressCalculator that implements stress calculations by accumulating each
timestep of stress using a decay equation. The decay aspect increases importance placed on
new damage. The decay equation is as follows:

`stress_accum = stress_accum(1-λ) + stress_instant(λ)

## Fields:
- `λ::AbstractFloat`: decay parameter used when calculating accumulated stress. Should be between 0 - 1.
"""
@kwdef struct DecayAreaScaledStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    λ::FT = 0.2

    function DecayAreaScaledStressCalculator{FT}(λ) where FT
        if λ < 0 || λ > 1
            @warn DECAY_ARG_WARNING
            λ = FT(0.2)
        end
        return new{FT}(λ)
    end
end

"""
    DecayAreaScaledStressCalculator(::Type{FT}; kwargs...)

A float type `FT` can be provided as the first argument of any
`DecayAreaScaledStressCalculator` constructor. A `DecayAreaScaledStressCalculator` of type
`FT` will be created by passing all other arguments to the correct constructor. 
"""
DecayAreaScaledStressCalculator(
    ::Type{FT},
    args...;
    kwargs...,
) where {FT <: AbstractFloat} =
    DecayAreaScaledStressCalculator{FT}(args...; kwargs...)

"""
    DecayAreaScaledStressCalculator(; kwargs...)

If type isn't specified as the first argument, `DecayAreaScaledStressCalculator`(; kwargs...) will
be of type `Float64` and the correct constructor will be called with all other arguments.
"""
DecayAreaScaledStressCalculator(args...; kwargs...) = DecayAreaScaledStressCalculator{Float64}(args...; kwargs...)

#= Updating `floe`'s accumulated and instantanious stress using a decay equation given the
current stress at a timestep `curr_stress::Matrix{FT}`.=#
function _update_stress!(stress_calculator::DecayAreaScaledStressCalculator, curr_stress, floe)
    λ = stress_calculator.λ
    floe.stress_accum .= (1 - λ) * floe.stress_accum + λ * curr_stress
    floe.stress_instant .= curr_stress
    return
end

"""
    DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}

Type of AbstractStressCalculator that calculates stress with each timestep with
`damage * stress_instant`, as suggested by Mukund Gupta. This method could keep track of
damage directly within each floe using the `damage` floe field, perhapes as a value between
0-1, and rather than calculating an "accumulated stress" to store damage, as done in
`DecayAreaScaledStressCalculator`.

In this calculator, the floe's `damage` field keeps track of an explicit parameter.
## Fields:

## Note:
This method is not implemented and thus throws an error upon creation. The infrastructure
for the damage calculator is provided, but functions that depend on this calculator need to
be implemented as detailed in the `AbstractStressCalculator` API. 
"""
@kwdef struct DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    function DamageStressCalculator{FT}() where FT
        throw(ErrorException("DamageStressCalculator isn't implemented yet!! Take a look at the issue/comments if you're interested in implementing."))
        return new{FT}()
    end
end

"""
    DamageStressCalculator(::Type{FT}; kwargs...)

A float type `FT` can be provided as the first argument of any
`DamageStressCalculator` constructor. A `DamageStressCalculator` of type `FT`
will be created by passing all other arguments to the correct constructor. 
"""
DamageStressCalculator(
    ::Type{FT},
    args...;
    kwargs...,
) where {FT <: AbstractFloat} =
    DamageStressCalculator{FT}(args...; kwargs...)

"""
    DamageStressCalculator(; kwargs...)

If type isn't specified, `DamageStressCalculator`(; kwargs...) will be of type
`Float64` and the correct constructor will be called with all other arguments.
"""
DamageStressCalculator(args...; kwargs...) = DamageStressCalculator{Float64}(args...; kwargs...)

#= Updating `floe`'s damage parameter depending on the current floe's instantanious stress
and then update the floe's accumulated stress based on this damage. =#
function _update_stress!(::DamageStressCalculator, curr_stress, floe)
    #= TODO: add in the changes to floe.damage and add any needed parameters to the
    DamageStressCalculator struct!! 
    floe.damage = ... =#
    floe.stress_accum .= floe.damage * curr_stress
    floe.stress_instant .= curr_stress
    return
end