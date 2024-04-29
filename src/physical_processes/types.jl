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