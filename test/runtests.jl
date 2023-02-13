using Subzero, LibGEOS, JLD2, NCDatasets, StructArrays, SplitApplyCombine, Statistics, VoronoiCells, GeometryBasics, Random, PolygonInbounds
using Test

@testset "Subzero.jl" begin
    #include("test_collisions.jl")
    #include("test_coupling.jl")
    include("test_floe.jl")
    include("test_floe_utils.jl")
    #include("test_model.jl")
    #include("test_output.jl")
end
