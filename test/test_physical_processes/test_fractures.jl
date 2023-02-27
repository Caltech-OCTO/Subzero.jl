@testset "Fractures" begin
    @testset "Fracture Criteria" begin
        # Test NoFracturee criteria
        @test NoFracture() isa NoFracture
        # Test HiblerYieldCurve criteria
        @test_throws ArgumentError HiblerYieldCurve(2.25, 20.0, [[[0.0, 0.0]]])
        @test HiblerYieldCurve(
            2.25,
            20.0,
            [[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]]
        ) isa HiblerYieldCurve
        # Test calculate_hibler
        hibler_verts = Subzero.calculate_hibler(0.5, 5e5, -1)
        hibler_poly = LibGEOS.Polygon(hibler_verts)
        @test isapprox(LibGEOS.area(hibler_poly), 49054437859.374, atol = -1e3)
        @test all(isapprox.(
            Subzero.find_poly_centroid(hibler_poly),
            [-1.25e5, -1.25e5],
            atol = 1e-3
        ))
        x_verts, y_verts = Subzero.seperate_xy(hibler_verts)
        @test all(isapprox.(
            extrema(x_verts),
            [-264743.588, 14727.999],
            atol = 1e-3
        ))
        @test all(isapprox.(
            extrema(y_verts),
            [-264743.588, 14727.999],
            atol = 1e-3
        ))
        hibler_verts = Subzero.calculate_hibler(0.25, 2.25e5, 20.0)
        hibler_poly = LibGEOS.Polygon(hibler_verts)
        @test isapprox(LibGEOS.area(hibler_poly), 2483380916.630, atol = -1e3)
        @test all(isapprox.(
            Subzero.find_poly_centroid(hibler_poly),
            [-28125, -28125],
            atol = 1e-3
        ))
        x_verts, y_verts = Subzero.seperate_xy(hibler_verts)
        @test all(isapprox.(
            extrema(x_verts),
            [-59567.307, 3313.799],
            atol = 1e-3
        ))
        @test all(isapprox.(
            extrema(y_verts),
            [-59567.307, 3313.799],
            atol = 1e-3
        ))
        # Test update criteria
        floes = StructArray([Floe(
            [[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]],
            0.25,  # Floe has a height of 0.25
            0.0,
        )])
        yield_curve = HiblerYieldCurve(floes)
        verts = deepcopy(yield_curve.vertices)
        @test yield_curve isa HiblerYieldCurve
        @test yield_curve.pstar == 2.25e5 && yield_curve.c == 20
        floes.height .= 0.5
        Subzero.update_criteria!(yield_curve, floes)
        @test verts != yield_curve.vertices

    end
    @testset "Fractures Floes" begin
        # Test determine_fractures

        # Test deform_floe!

        # Test split_floe

        # Test fracture_floes!

    end
end