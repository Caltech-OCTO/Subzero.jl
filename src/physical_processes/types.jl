"""
Structs and functions for fracturing floes
"""

"""
    AbstractStressCalculator

Abstract type ways of keeping track of stress.
"""
abstract type AbstractStressCalculator end

"""
    RunningAverageCalculator

idk yet.
"""
@kwdef struct RunningAverageCalculator <: AbstractStressCalculator
    nhistory::Int = 100

    function RunningAverageCalculator(nhistory) 
        if nhistory < 1
            @warn "Value of nhistory must be greater than or equal to 1. Resetting to default value of 100."
            nhistory = 100
        end

        return new(nhistory)
    end
end

# TODO: Get rid of this calculator
"""
    DecayCalculator

idk yet.
"""
@kwdef struct DecayCalculator{FT<:AbstractFloat} <: AbstractStressCalculator
    τ::FT = 20.0

    function DecayCalculator{FT}(τ::FT) where {FT<:AbstractFloat}
        # if τ > 1
        #     @warn "Value of τ must be less than or equal to 1. Resetting to default value of 0.1."
        #     τ = 20.0
        # end

        return new(τ)
    end
end

DecayCalculator(args...; kwargs...) = DecayCalculator{Float64}(args...; kwargs...)

# Rename this DecayAreaScaledCalc...
"""
    AreaScaledCalculator

idk yet.
"""
@kwdef struct AreaScaledCalculator{FT<:AbstractFloat} <: AbstractStressCalculator
    τ::FT = 20.0
    α::FT = 0.5

    function AreaScaledCalculator{FT}(τ::FT, α::FT) where {FT<:AbstractFloat}
        # if τ > 1
        #     @warn "Value of τ must be less than or equal to 1. Resetting to default value of 0.1."
        #     τ = 20.0
        # end

        return new{FT}(τ, α)
    end
end

AreaScaledCalculator(args...; kwargs...) = AreaScaledCalculator{Float64}(args...; kwargs...)

# TODO: New Calculator: DamageStressCalc (Mukund's)
"""
    DamageStressCalculator

idk yet.
"""
@kwdef struct DamageStressCalculator{FT<:AbstractFloat} <: AbstractStressCalculator
    τ::FT = 20.0

    function DamageStressCalculator{FT}(τ::FT) where {FT<:AbstractFloat}
        # if τ > 1
        #     @warn "Value of τ must be less than or equal to 1. Resetting to default value of 0.1."
        #     τ = 20.0
        # end

        return new(τ)
    end
end

DamageStressCalculator(args...; kwargs...) = DamageStressCalculator{Float64}(args...; kwargs...)