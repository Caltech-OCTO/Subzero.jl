@testset "Update floe" begin
    @testset "Stress/Strain" begin
        # Test Stress History buffer
        scb = Subzero.StressCircularBuffer{Float64}(10)
        @test scb.cb.capacity == 10
        @test scb.cb.length == 0
        @test eltype(scb.cb) <: Matrix{Float64}
        @test scb.total isa Matrix{Float64}
        @test all(scb.total .== 0)
        fill!(scb, ones(Float64, 2, 2))
        @test scb.cb.length == 10
        @test all(scb.total .== 10)
        @test scb.cb[1] == scb.cb[end] == ones(Float64, 2, 2)
        @test mean(scb) == ones(Float64, 2, 2)
        push!(scb, zeros(Float64, 2, 2))
        @test scb.cb.length == 10
        @test all(scb.total .== 9)
        @test scb.cb[end] == zeros(Float64, 2, 2)
        @test scb.cb[1] == ones(Float64, 2, 2)
        @test all(mean(scb) .== 0.9)

        # Test Stress and Strain Calculations
        floe_dict = load(
            "inputs/test_values.jld2"  # uses the first 2 element
        )
        stresses = [[-10.065, 36.171, 36.171, -117.458],
            [7.905, 21.913, 21.913, -422.242]]
        stress_histories = [[-4971.252, 17483.052, 17483.052, -57097.458],
            [4028.520, 9502.886, 9502.886, -205199.791]]
        strains = [[-3.724, 0, 0, 0], [7.419, 0, 0,	-6.987]]
        strain_multiplier = [1e28, 1e6]

        for i in 1:2
            f = Floe(
                floe_dict["coords"][i],
                floe_dict["height"][i],
                0.0,
                u = floe_dict["u"][i],
                v = floe_dict["v"][i],
                ξ = floe_dict["ξ"][i],
            )
            f.interactions = floe_dict["interactions"][i]
            f.stress_history = floe_dict["stress_history"][i]
            stress = Subzero.calc_stress!(f)
            @test all(isapprox.(vec(f.stress), stresses[i], atol = 1e-3))
            @test all(isapprox.(
                vec(f.stress_history.cb[end]),
                stress_histories[i],
                atol = 1e-3
            ))
            Subzero.calc_strain!(f)
            @test all(isapprox.(
                vec(f.strain) .* strain_multiplier[i],
                strains[i],
                atol = 1e-3
            ))
        end
    end
    @testset "Floe shape" begin
        # Test replace floe
        coords1 = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [10.0, 0.0],
            [0.0, 0.0],
        ]]
        f1 = Floe(coords1, 0.5, 0.0)  # this is a square
        mass1 = f1.mass
        tri_coords = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [0.0, 0.0],
        ]]
        tri_poly = LG.Polygon(tri_coords)  # this is a triangle
        Subzero.replace_floe!(
            f1,
            tri_poly,
            f1.mass,
            Constants(),
            100,
            Xoshiro(1)
        )
        @test f1.centroid == Subzero.find_poly_centroid(tri_poly)
        @test f1.coords == Subzero.find_poly_coords(tri_poly)
        @test f1.area == LG.area(tri_poly)
        @test f1.mass == mass1
        @test f1.height * f1.area * 920.0 == f1.mass
        @test f1.α == 0
        @test f1.status.tag == Subzero.active
        @test f1.rmax == 10*√5 / 3
    end
end