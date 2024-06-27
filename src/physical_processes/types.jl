"""
Structs and functions for calculating stress
"""

"""
    AbstractStressCalculator

Abstract type ways of keeping track of stress.
"""
abstract type AbstractStressCalculator{FT<:AbstractFloat}  end

"""
    DecayAreaScaledCalculator{FT<:AbstractFloat}<:AbstractStressCalculator

Type of AbstractStressCalculator that implements stress calculations by accumulating each
timestep of stress using a decay equation. The decay aspect increases importance placed on
new damage. The decay equation is as follows:

stress_accum = stress_accum + Δt/τ .* (stress_instant - stress_accum)

To avoid recomputation, Δt/τ is saved in the floe's damage field. 

The area-scaling part of the stress calculations comes into play when deciding if a floe
will be fractured. The accumulated stress can be scaled by a ratio of
(floe.area/min_floe_area)^α. By changing the α value, either larger or smaller floes
fracture more eassily. The scaled value will not be saved as this is equivalent to scaling
the simulation fracture criteria (morh's cone or hibler's ellipse, etc.), but it is less
computationally intensive. 

Fields:
    λ      <AbstractFloat> Decay parameter used when calculating accumulated stress in
                update_stress! function. Commonly set to Δt/τ. Should be between 0 - 1 such
                that stress_accum = stress_accum(1-λ) + stress_instant(λ)
    α      <AbstractFloat> Adjusts ellipse in stress space by raising area ratio to the α
Note:
    λ is used in calc_stress!, whereas α is used in determine_fractures.
"""
@kwdef struct DecayAreaScaledCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    λ::FT = 0.2
    α::FT = 0.0
end

"""
    DecayAreaScaledCalculator(::Type{FT}; kwargs...)

A float type FT can be provided as the first argument of any
DecayAreaScaledCalculator constructor. A DecayAreaScaledCalculator of type FT
will be created by passing all other arguments to the correct constructor. 
"""
DecayAreaScaledCalculator(
    ::Type{FT},
    args...;
    kwargs...,
) where {FT <: AbstractFloat} =
    DecayAreaScaledCalculator{FT}(args...; kwargs...)

"""
    DecayAreaScaledCalculator(; kwargs...)

If type isn't specified, DecayAreaScaledCalculator(; kwargs...) will be of type
Float64 and the correct constructor will be called with all other arguments.
"""
DecayAreaScaledCalculator(args...; kwargs...) = DecayAreaScaledCalculator{Float64}(args...; kwargs...)

"""
    DamageStressCalculator

Type of AbstractStressCalculator that calculates stress with damage * stress_instant, as
suggested by Mukund Gupta. This method could keep track of damage directly, perhapes as a
value between 0-1, and rather than calculating an "accumulated stress" to store damage, as
done in DecayAreaScaledCalculator.

In this calculator, the floe's damage field keeps track of an explicit parameter for damage.
Fields:
    τ      <AbstractFloat> Difference between current and previous stress scaled by Δt/τ.
                           This field should represent damage.
Note:
    This method is not fully implemented. The infrastructure for the damage calculator is
    provided, but functions that depend on this calculator need to be implemented. 
    There is no α parameter in this calculator because currently it does not adjust the
    boundary in stress space by multiplying the eigenvalues of stress_accum by something.
    This could be implemented if the user desires.

    The functions that need to be implemented are: scale_stress!, update_damage!, and
    update_stress!. You may also need to create an inialize damage parameter, which can set
    an initial damage based on the calculator type if you don't want it to start at 0.

    You can add needed fields to the below struct as in the DecayAreaScaledCalculator
"""
@kwdef struct DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    function DamageStressCalculator{FT}() where FT
        throw(ErrorException("DamageStressCalculator isn't implemented yet!! Take a look at the issue/comments if you're interested in implementing."))
        return new{FT}()
    end

end

"""
    DamageStressCalculator(::Type{FT}; kwargs...)

A float type FT can be provided as the first argument of any
DamageStressCalculator constructor. A DamageStressCalculator of type FT
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

If type isn't specified, DamageStressCalculator(; kwargs...) will be of type
Float64 and the correct constructor will be called with all other arguments.
"""
DamageStressCalculator(args...; kwargs...) = DamageStressCalculator{Float64}(args...; kwargs...)