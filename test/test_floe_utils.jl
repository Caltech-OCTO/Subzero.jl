@testset "Floe Utils" begin
    ext = [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    hole1 = [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3], [0.2, 0.3]]
    hole2 = [[0.5, 0.6], [0.5, 0.5], [0.6, 0.5], [0.6, 0.6], [0.5, 0.6]]
    poly_nohole = Subzero.make_polygon([ext])
    poly_hole1 = Subzero.make_polygon([ext, hole1])
    poly_hole2 = Subzero.make_polygon([ext, hole1, hole2])
    
    # Test validating/correcting RingVecs and PolyVecs
    @test Subzero.valid_ringvec!(ext) == ext
    invalid_ext = [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
    @test Subzero.valid_ringvec!(invalid_ext) == ext
    @test_throws AssertionError Subzero.valid_ringvec!([[0.0, 1.0], [0.0, 0.0]])
    invalid_coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
                      [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3]]]
    @test Subzero.valid_polyvec!(invalid_coords) == [ext, hole1]
    duplicate_invalid_ext = [
        [0.0, 1.0],
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 0.0], 
        [1.0, 1.0],
    ]
    @test Subzero.valid_ringvec!(duplicate_invalid_ext) == ext
    @test Subzero.valid_ringvec!(
        [[0.0, 1.0], [0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
    ) == ext
    @test Subzero.valid_ringvec!(
        [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
    ) == ext
    @test_throws AssertionError Subzero.valid_ringvec!(
        [[0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
    )
    @test_throws AssertionError Subzero.valid_polyvec!([[Float64[]]])

    # Test predicate hashole for coords, polygons and multipolygons
    @test !Subzero.hashole([ext])
    @test Subzero.hashole([ext, hole1])
    @test !Subzero.hashole(poly_nohole)
    @test Subzero.hashole(poly_hole1)
    @test Subzero.hashole(poly_hole2)

    # Test removing holes from polygons
    copy_holes = [ext, hole1]
    poly_copy_holes = Subzero.make_polygon(copy_holes)
    Subzero.rmholes!(copy_holes)
    @test copy_holes == [ext]
    Subzero.rmholes!(poly_copy_holes)
    @test GI.nhole(poly_copy_holes) == 0

    # Test translating coordinates and polygons
    @test Subzero.translate([ext], 0.0, 0.0) == [ext]
    trans_ext = Subzero.translate([ext], 1.0, 2.0)
    @test trans_ext == [[[1.0, 3.0],  [1.0, 2.0],  [2.0, 2.0],
                        [2.0, 3.0], [1.0, 3.0]]]
    copy_ext = [deepcopy(ext)]
    Subzero.translate!(copy_ext, 1.0, 2.0)
    @test copy_ext == trans_ext
    test_trans = [[[-2.0, 2.0], [-2.0, 1.0], [-1.0, 1.0], [-1.0, 2.0]]]
    @test Subzero.translate(test_trans, 1.5, -1.5) ==
        [[[-0.5, 0.5], [-0.5, -0.5], [0.5, -0.5], [0.5, 0.5]]]
    Subzero.translate!(test_trans, 1.5, -1.5)
    @test test_trans == [[[-0.5, 0.5], [-0.5, -0.5], [0.5, -0.5], [0.5, 0.5]]]

    # Test moment of intertia calculations - compared to values output my MATLAB
    poly_moment = Subzero._calc_moment_inertia(Float64, Subzero.make_polygon([ext]), [0.5, 0.5], 0.25)
    @test isapprox(poly_moment, 38.333, atol = 0.001)
    @test Subzero._calc_moment_inertia(Float64, poly_nohole, GO.centroid(poly_nohole), 0.25) == poly_moment
    tri_poly = Subzero.make_polygon([[[0, 1], [0, 0], [1, 0], [0, 1]]] .* 6.67)
    tri_moment = Subzero._calc_moment_inertia(Float64, tri_poly, GO.centroid(tri_poly), 0.5)
    @test isapprox(tri_moment, 50581.145, atol = 0.001)

    # ------------------ Test rotating coordinates ------------------
    og_coords = [
        [[-1.0, -1.0], [-1, 1], [1, 1], [1, -1], [-1, -1]],
        [[-1.0, -1.0], [-1, 1], [1, 1], [1, -1], [-1, -1]],
    ]
    copy_coords = deepcopy(og_coords)
    Subzero.rotate_radians!(copy_coords, π/4)
    same = true
    answer = [[[0, -√2], [-√2, 0], [0, √2], [√2, 0], [0, -√2]],
        [[0, -√2], [-√2, 0], [0, √2], [√2, 0], [0, -√2]]]
    for i in eachindex(copy_coords)
        same = same && all(isapprox.(copy_coords[i], answer[i]))
    end
    @test same
    Subzero.rotate_radians!(copy_coords, 7π/4)
    same = true
    for i in eachindex(copy_coords)
        same = same & all(isapprox.(copy_coords[i], og_coords[i]))
    end
    @test same
end