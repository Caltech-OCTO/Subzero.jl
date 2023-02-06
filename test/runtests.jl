using Subzero, LibGEOS, JLD2, NCDatasets, StructArrays, SplitApplyCombine, Statistics
using Test

@testset "Subzero.jl" begin
    include("test_floe_operations.jl")
    include("test_model.jl")
    include("test_coupling.jl")
    include("test_collisions.jl")
    #include("test_output.jl")
    #include("test_simulation.jl")
end
