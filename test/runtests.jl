using DataStructures, Subzero, JLD2, NCDatasets, StructArrays,
    SplitApplyCombine, Statistics, VoronoiCells, GeometryBasics, Random,
    PolygonInbounds
import LibGEOS as LG
using Test

@testset "Subzero.jl" begin
    include("test_physical_processes/test_update_floe.jl")
    include("test_physical_processes/test_collisions.jl")
    include("test_physical_processes/test_coupling.jl")
    include("test_physical_processes/test_process_settings.jl")
    include("test_physical_processes/test_fractures.jl")
    include("test_physical_processes/test_simplification.jl")
    include("test_physical_processes/test_ridge_raft.jl")
    include("test_floe.jl")
    include("test_floe_utils.jl")
    include("test_model.jl")
    include("test_output.jl")
    include("test_simulation.jl")
    include("test_conservation.jl")
end
