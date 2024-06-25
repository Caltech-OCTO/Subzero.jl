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
    MovingBoundary,
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
    MohrsCone,
    CollisionSettings,
    FractureSettings,
    CouplingSettings,
    SimplificationSettings,
    RidgeRaftSettings,
    WeldSettings,
    FloeSettings,
    PolyVec,
    OutputWriters,
    check_energy_momentum_conservation_julia,
    IceStressCell,
    CellFloes,
    MonteCarloPointsGenerator,
    SubGridPointsGenerator,
    plot_sim

import Base.@kwdef # this is being exported as of version 1.9
import GeometryOps as GO
import GeometryOps.GeoInterface as GI
import StaticArrays as SA
using CairoMakie, CoordinateTransformations, DataStructures, Dates, Extents,
    Interpolations, JLD2, LinearAlgebra, Logging, Makie, Measures, NCDatasets, NetCDF,
    Printf, Random, Rotations, SplitApplyCombine, Statistics, StructArrays,
    VoronoiCells

"""
Coordinates are vector of vector of vector of points of the form:
[[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]], 
 [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]], ...] where the xy coordinates
 are the exterior border of the floe and the wz coordinates, or any other
 following sets of coordinates, describe holes within the floe.
 This form is for easy conversion to polygons.
"""
const PolyVec{T} = Vector{Vector{Vector{T}}} where T<:Real

"""
Coordinates are vector of vector of points of the form:
[[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]] where the xy coordinates form a
closed ring. PolyVec objects can be made out RingVec objects.
This form is for each conversion to LinearRings, which can also be made into Polygons.
"""
const RingVec{T} = R where {
    T<:Real,
    V<:AbstractArray{T},
    R <: AbstractArray{V},
}

const Polys{T} = GI.Polygon{false, false, Vector{GI.LinearRing{false, false, Vector{Tuple{T, T}}, Nothing, Nothing}}, Nothing, Nothing} where T <: AbstractFloat

Base.convert(::Type{Polys{Float32}}, p::Polys{Float64}) = GO.tuples(p, Float32)
Base.convert(::Type{Polys{Float64}}, p::Polys{Float32}) = GO.tuples(p, Float64)

# Model
include("simulation_components/floe.jl")
include("floe_utils.jl")
include("simulation_components/domain_and_grid.jl")
include("simulation_components/ocean_and_atmos.jl")
include("simulation_components/model.jl")
# Physical Processes
include("physical_processes/update_floe.jl")
include("physical_processes/coupling.jl")
include("physical_processes/collisions.jl")
include("physical_processes/fractures.jl")
include("physical_processes/process_settings.jl")
include("physical_processes/simplification.jl")
include("physical_processes/ridge_raft.jl")
include("physical_processes/welding.jl")
# Tools
include("tools/plotting.jl")
include("tools/conservation_em.jl")
include("tools/compare_files.jl")
include("logger.jl")
# Simulation
include("output.jl")
include("simulation_components/simulation.jl")
end