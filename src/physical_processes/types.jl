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

"""
    DecayCalculator

idk yet.
"""
@kwdef struct DecayCalculator <: AbstractStressCalculator
    τ::Float32 = 0.1

    function DecayCalculator( τ) 
        if τ > 1
            @warn "Value of τ must be less than or equal to 1. Resetting to default value of 0.1."
            τ = 0.1
        end

        return new(τ)
    end
end