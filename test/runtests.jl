using Subzero, LibGEOS, JLD2, StructArrays, SplitApplyCombine
using Test

@testset "Subzero.jl" begin
    #include("test_floe_operations.jl")
    #include("test_model.jl")
    #include("test_coupling.jl")
    include("test_collisions.jl")
    #include("test_simulation.jl")
end
