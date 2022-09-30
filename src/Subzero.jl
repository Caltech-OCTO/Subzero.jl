"""
Module `Subzero.jl` - UW's sea ice model ported from MATLAB to Julia
"""
module Subzero

export
    Grid,
    Ocean,
    Wind,
    Floe,
    Model,
    Simulation,
    run!

import LibGEOS as LG
import Base.@kwdef # this is being exported as of version 1.9
using Plots, StructArrays, Statistics


const PolyVec{T} = Vector{Vector{Vector{T}}} where T<:AbstractFloat


# some structs here!

include("floe_operations.jl")
include("model.jl")

end
