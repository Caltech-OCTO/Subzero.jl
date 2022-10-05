using Subzero, LibGEOS
using Test

@testset "Subzero.jl" begin
    include("test_floe_operations.jl")
    include("test_model.jl")
end
