@testset "Simulation" begin
    @testset "Stress/Strain" begin
        floes = load("inputs/test_floes.jld2", "stress_strain_floe1", "stress_strain_floe2")
        stresses = [[-10.065, 36.171, 36.171, -117.458], [7.905, 21.913, 21.913, -422.242]]
        stress_histories = [[-4971.252, 17483.052, 17483.052, -57097.458], [4028.520, 9502.886, 9502.886, -205199.791]]
        strains = [[-3.724, 0, 0, 0], [7.419, 0, 0,	-6.987]]
        strain_multiplier = [1e28, 1e6]

        for i in eachindex(floes)
            f = floes[i]
            Subzero.calc_stress_strain!(f)
            @test all(isapprox.(vec(f.stress), stresses[i], atol = 1e-3))
            @test all(isapprox.(vec(f.stress_history[end]), stress_histories[i], atol = 1e-3))
            @test all(isapprox.(vec(f.strain) .* strain_multiplier[i], strains[i], atol = 1e-3))
        end
    end
end