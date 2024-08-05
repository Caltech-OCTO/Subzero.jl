using Test

@testset "Fracture Criteria" begin
    FT = Float64
    mean_height = 0.5
    pstar, c, compactness = 5e5, -1.0, 1.0
    # Test NoFracturee criteria
    @test NoFracture() isa NoFracture
    # Test HiblerCurveFractureCriteria criteria
    # @test HiblerCurveFractureCriteria(;
    #     pstar = 2.25,
    #     c = 20.0,
    #     poly = Subzero.make_polygon([[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]]),
    # ) isa HiblerCurveFractureCriteria
    # Test _calculate_hibler
    hibler_a, hibler_b, hibler_h = Subzero._calculate_hibler(FT, mean_height, pstar, c, compactness)
    # @test isapprox(GO.area(hibler_poly), 49054437859.374, atol = -1e3)
    # @test all(isapprox.(
    #     GO.centroid(hibler_poly),
    #     (-1.25e5, -1.25e5),
    #     atol = 1e-3
    # ))
    # hibler_verts = Subzero.find_poly_coords(hibler_poly)
    # x_verts, y_verts = first.(hibler_verts[1]), last.(hibler_verts[1])
    # @test all(isapprox.(
    #     extrema(x_verts),
    #     [-264743.588, 14727.999],
    #     atol = 1e-3
    # ))
    # @test all(isapprox.(
    #     extrema(y_verts),
    #     [-264743.588, 14727.999],
    #     atol = 1e-3
    # ))
    hibler_a, hibler_b, hibler_h = Subzero._calculate_hibler(FT, 0.25, 2.25e5, 20.0, 1.0)
    # hibler_verts = Subzero.find_poly_coords(hibler_poly)
    # @test isapprox(GO.area(hibler_poly), 2483380916.630, atol = -1e3)
    # @test all(isapprox.(
    #     GO.centroid(hibler_poly),
    #     (-28125, -28125),
    #     atol = 1e-3
    # ))
    # x_verts, y_verts = first.(hibler_verts[1]), last.(hibler_verts[1])
    # @test all(isapprox.(
    #     extrema(x_verts),
    #     [-59567.307, 3313.799],
    #     atol = 1e-3
    # ))
    # @test all(isapprox.(
    #     extrema(y_verts),
    #     [-59567.307, 3313.799],
    #     atol = 1e-3
    # ))
    # Test Mohr's Cone
    @test typeof(Subzero.MohrsConeFractureCriteria(Float32)) <: MohrsConeFractureCriteria{Float32}
    @test typeof(Subzero.MohrsConeFractureCriteria(Float64)) <: MohrsConeFractureCriteria{Float64}

    # Float64 Mohr's Cone with q, σc, σ11
    mohrs_verts_64 = Subzero.find_poly_coords(Subzero._calculate_mohrs(FT, 5.2, 2.5e5, -3.375e4))
    @test all(isapprox.(
        mohrs_verts_64[1],
        [
            [59523.809, 59523.809],
            [33750.0, -74500.0],
            [-74500.0, 33750.0],
            [59523.809, 59523.809],
        ],
        atol = 1e-3
    ))
    # Float32 Mohr's Cone with q, σc, σ11
    mohrs_verts_32 = Subzero.find_poly_coords(Subzero._calculate_mohrs(FT, 5.2, 2.5e5, 1.5e5))
    @test all(isapprox.(
        mohrs_verts_32[1],
        [
            [59523.809, 59523.809],
            [-150000.0, -1.03e6],
            [-1.03e6, -150000.0],
            [59523.809, 59523.809],
        ],
        atol = 1e-3,
    ))
    # Float64 Mohr's Cone with σ1, σ2, σ11, σ22
    mohrs_verts_coords = Subzero._calculate_mohrs(
        FT,
        5.95e4,
        5.95e4,
        -1.5e5,
        -1.03e6,
    )
    # Test update criteria for Hibler's curve
    floes = StructArray([Floe(
        [[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]],
        0.25,  # Floe has a height of 0.25
        0.0,
    )])
    yield_curve = HiblerCurveFractureCriteria(; floes)
    # old_poly = yield_curve.poly
    @test yield_curve isa HiblerCurveFractureCriteria
    @test yield_curve.pstar == 2.25e5 && yield_curve.c == 20
    floes.height .= 0.5
    Subzero._update_criteria!(yield_curve, floes)
    # @test !GO.equals(old_poly, yield_curve.poly)
    # Test update criteria for Mohr's cone
    cone_curve = MohrsConeFractureCriteria()
    # old_poly = cone_curve.poly
    Subzero._update_criteria!(cone_curve, floes)
    # @test GO.equals(old_poly, cone_curve.poly)
end