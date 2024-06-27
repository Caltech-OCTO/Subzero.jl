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

Type of AbstractStressCalculator that implements stress calculations the same way that 
Brandon Montemuro and Georgy Manucharyan do in the MatLab versio of the code. The area-scaling
part of the stress calculations keep a record of previous stress (and implicity damage). The
decay aspect increases importance placed on new damage.
Fields:
    τ      <AbstractFloat> Difference between current and previous stress scaled by Δt/τ
    α      <AbstractFloat> Adjusts ellipse in stress space by raising area ratio to the α
Note:
    τ is used in calc_stress!(), whereas α is used in determine_fractures().
"""
@kwdef struct DecayAreaScaledCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    τ::FT = 25.0
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

Type of AbstractStressCalculator that calculates stress with damage*currstress, as suggested
by Mukund Gupta. This method keeps track of damage directly and rather than modifying the 
"accumulated stress" to store damage, as done in DecayAreaScaledCalculator. This calculator
keeps track of an explicit parameter for damage.
Fields:
    τ      <AbstractFloat> Difference between current and previous stress scaled by Δt/τ.
                           This field should represent damage.
Note:
    This method is not fully implemented. The infrastructure for the damage calculator is
    provided, but functions that depend on this calculator,  
    There is no α parameter in this calculator because currently it does not adjust the
    boundary in stress space by multiplying the eigenvalues of stress_accum by something.
    This could be implemented if the user desires.
"""
@kwdef struct DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator{FT}
    τ::FT = 20.0
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