"""
Module `Subzero.jl` - UW's sea ice model ported from MATLAB to Julia
"""
module Subzero

export
    AbstractGrid,
    RegRectilinearGrid,
    Ocean,
    Atmos,
    OpenBC,
    PeriodicBC,
    CollisionBC,
    CompressionBC,
    AbstractBoundary,
    CompressionBoundary,
    CollisionBoundary,
    PeriodicBoundary,
    OpenBoundary,
    North,
    South,
    East,
    West,
    Domain,
    TopographyElement,
    Floe,
    Constants,
    Model,
    Simulation,
    timestep_sim!,
    run!,
    AbstractOutputWriter,
    CheckpointOutputWriter,
    GridOutputWriter, 
    FloeOutputWriter,
    InitialStateOutputWriter,
    GridOutput,
    FloeOutput,
    add_ghosts!,
    # Interaction field enum and elements
    InteractionFields,
    floeidx,
    xforce,
    yforce,
    xpoint,
    ypoint,
    torque,
    overlap, 
    initialize_floe_field,
    initialize_topography_field,
    AbstractFractureCriteria,
    NoFracture,
    HiblerYieldCurve,
    CollisionSettings,
    FractureSettings,
    CouplingSettings,
    SimplificationSettings,
    PolyVec,
    OutputWriters,
    check_energy_momentum_conservation_julia,
    IceStressCell,
    CellFloes

import Base.@kwdef # this is being exported as of version 1.9
import LibGEOS as LG
using DataStructures, GeometryBasics, Interpolations, JLD2, LinearAlgebra,
    Logging, Measures, NamedArrays, NCDatasets, NetCDF, Plots, PolygonInbounds,
    Printf, Random, SplitApplyCombine, StaticArrays, Statistics, StructArrays,
    VoronoiCells

"""
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates
 are the exterior border of the floe and the wz coordinates, or any other
 following sets of coordinates, describe holes within the floe.
 This form is for easy conversion to LibGEOS Polygons.
"""
const PolyVec{T} = Vector{Vector{Vector{T}}} where T<:Real

"""
Coordinates are vector of vector of points of the form:
[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]] where the xy coordinates form a
closed ring. PolyVec objects can be made out RingVec objects.
This form is for each conversion to LibGEOS LinearRings, which can also be made into Polygons.
"""
const RingVec{T} = R where {
    T<:Real,
    V<:AbstractArray{T},
    R <: AbstractArray{V},
}

# Model
include("floe.jl")
include("floe_utils.jl")
include("model.jl")
# Physical Processes
include("physical_processes/update_floe.jl")
include("physical_processes/coupling.jl")
include("physical_processes/collisions.jl")
include("physical_processes/fractures.jl")
include("physical_processes/process_settings.jl")
include("physical_processes/simplification.jl")
# Tools
include("tools/plotting.jl")
include("tools/conservation_em.jl")
include("tools/compare_files.jl")
# Simulation
include("output.jl")
include("simulation.jl")
end