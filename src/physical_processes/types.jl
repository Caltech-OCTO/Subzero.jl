"""
Structs and functions for fracturing floes
"""

"""
    AbstractStressCalculator

Abstract type ways of keeping track of stress.
"""
abstract type AbstractStressCalculator end

"""
    DecayAreaScaledCalculator{FT<:AbstractFloat}<:AbstractStressCalculator

Type of AbstractStressCalculator that implements stress calculations the same way that 
Brandon Montemuro and Georgy Manucharyan do in the MatLab version.
Fields:
    τ      <AbstractFloat> Difference between current and previous stress scaled by Δt/τ
    α      <AbstractFloat> Adjusts ellipse in stress space by raising area ratio to the α
Note:
    τ is used in calc_stress!(), whereas α is used in determine_fractures().
"""
@kwdef struct DecayAreaScaledCalculator{FT<:AbstractFloat} <: AbstractStressCalculator
    τ::FT = 20.0
    α::FT = 0.5

    function DecayAreaScaledCalculator{FT}(τ::FT, α::FT) where {FT<:AbstractFloat}
        # Can insert warnings about parameter values here
        return new{FT}(τ, α)
    end
end

DecayAreaScaledCalculator(args...; kwargs...) = DecayAreaScaledCalculator{Float64}(args...; kwargs...)

"""
    DamageStressCalculator

Type of AbstractStressCalculator that calculates stress with damage*curre_stress, as suggested
by Mukund Gupta.
    τ      <AbstractFloat> Difference between current and previous stress scaled by Δt/τ
Note:
    There is no α parameter in this calculator because currently it does not adjust the
    boundary in stress space by multiplying the eigenvalues of stress_accum by something.
    This could be implemented if the user desires.
"""
@kwdef struct DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator
    τ::FT = 20.0

    function DamageStressCalculator{FT}(τ::FT) where {FT<:AbstractFloat}
        # Can insert parameter warnings here

        return new(τ)
    end
end

DamageStressCalculator(args...; kwargs...) = DamageStressCalculator{Float64}(args...; kwargs...)