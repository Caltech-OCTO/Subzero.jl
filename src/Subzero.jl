"""
Module `Subzero.jl` - UW's sea ice model ported from MATLAB to Julia
"""
module Subzero

export
    Grid,
    Ocean,
    Wind,
    AbstractBC,
    Open,
    Boundary,
    AbstractDomain,
    RectangleDomain,
    Floe,
    Model,
    Simulation,
    run!

import LibGEOS as LG
import Base.@kwdef # this is being exported as of version 1.9
using Plots, StructArrays, Statistics


"""
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates are the exterior border of the floe and the wz coordinates, or any other following sets of coordinates, describe holes within the floe. This form is for easy conversion to LibGEOS Polygons.
"""
const PolyVec{T} = Vector{Vector{Vector{T}}} where T<:AbstractFloat
"""
Coordinates are vector of vector of points of the form:
[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]] where the xy coordinates form a closed ring. PolyVec objects can be made out RingVec objects. This form is for each conversion to LibGEOS LinearRings, which can also be made into Polygons.
"""
const RingVec{T} = Vector{Vector{T}} where T<:AbstractFloat

include("floe_operations.jl")
include("model.jl")

end
