@testset "Floe Operations" begin
    ext = [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    hole1 = [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3], [0.2, 0.3]]
    hole2 = [[0.5, 0.6], [0.5, 0.5], [0.6, 0.5], [0.6, 0.6], [0.5, 0.6]]
    poly_nohole = LibGEOS.Polygon([ext])
    poly_hole1 = LibGEOS.Polygon([ext, hole1])
    poly_hole2 = LibGEOS.Polygon([ext, hole1, hole2])
    multipoly_hole1 = LibGEOS.MultiPolygon([poly_nohole, poly_hole1])
    multipoly_hole2 = LibGEOS.MultiPolygon([poly_nohole, poly_hole1,
                                                         poly_hole2])

    # Test validating/correcting RingVecs and PolyVecs
    @test Subzero.valid_ringvec!(ext) == ext
    invalid_ext = [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
    @test Subzero.valid_ringvec!(invalid_ext) == ext
    @test_throws AssertionError Subzero.valid_ringvec!([[0.0, 1.0], [0.0, 0.0]])
    invalid_coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
                      [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3]]]
    @test Subzero.valid_polyvec!(invalid_coords) == [ext, hole1]

    # Test predicate hashole for coords, polygons and multipolygons
    @test !Subzero.hashole([ext])
    @test Subzero.hashole([ext, hole1])
    @test !Subzero.hashole(poly_nohole)
    @test Subzero.hashole(poly_hole1)
    @test Subzero.hashole(poly_hole2)
    @test Subzero.hashole(multipoly_hole1)

    # Test removing holes from polygons and multipolygons
    @test [ext] == Subzero.rmholes([ext])
    @test [ext] == Subzero.rmholes([ext, hole1])
    @test LibGEOS.equals(Subzero.rmholes(poly_nohole), poly_nohole)
    @test LibGEOS.equals(Subzero.rmholes(poly_hole1), poly_nohole)
    @test LibGEOS.equals(Subzero.rmholes(poly_hole2), poly_nohole)
    @test !Subzero.hashole(Subzero.rmholes(multipoly_hole1))

    # Test sorting regions from polygons and multipolygons
    @test Subzero.sortregions(poly_nohole) == [poly_nohole]
    poly_lst = Subzero.sortregions(multipoly_hole2)
    @test LibGEOS.equals(poly_lst[1], poly_nohole)
    @test LibGEOS.equals(poly_lst[3], poly_hole2)

    # Test translating coordinates and polygons
    @test Subzero.translate([ext], [0.0, 0.0]) == [ext]
    trans_ext = Subzero.translate([ext], [1.0, 2.0])
    @test trans_ext == [[[1.0, 3.0],  [1.0, 2.0],  [2.0, 2.0],
                        [2.0, 3.0], [1.0, 3.0]]]
    trans_nohole = Subzero.translate(poly_nohole, [1, 2])
    @test LibGEOS.equals(trans_nohole, LibGEOS.Polygon(trans_ext))
    trans_hole1 = Subzero.translate(poly_hole1, [1.0, 2.0])
    @test LibGEOS.equals(trans_hole1, LibGEOS.Polygon(
        [[[1.0, 3.0],  [1.0, 2.0],  [2.0, 2.0],  [2.0, 3.0], [1.0, 3.0]],
         [[1.2, 2.3], [1.2, 2.2], [1.3, 2.2], [1.3, 2.3], [1.2, 2.3]]]))

    # Test scaling polygons
    centered_coords = [[[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], 
                        [-1.0, 1.0], [-1.0, -1.0]]]
    centered_poly = LibGEOS.Polygon(centered_coords)
    @test LibGEOS.equals(Subzero.scale(centered_poly, 1), centered_poly)
    @test LibGEOS.equals(Subzero.scale(centered_poly, 2.0),
        LibGEOS.Polygon(centered_coords .* 2.0))

    # Test seperating x and y coordiantes from PolyVec form
    x_ext, y_ext = Subzero.seperate_xy([ext])
    @test x_ext == [0.0, 0.0, 1.0, 1.0, 0.0]
    @test y_ext == [1.0, 0.0, 0.0, 1.0, 1.0]

    # Test moment of intertia calculations - compared to values output my MATLAB
    poly_moment = Subzero.calc_moment_inertia([ext], 0.25)
    @test isapprox(poly_moment, 153.333, atol = 0.001)
    @test Subzero.calc_moment_inertia(poly_nohole, 0.25) == poly_moment
    @test Subzero.calc_moment_inertia(poly_hole1, 0.25) == poly_moment
    tri_coords = [[[0, 1], [0, 0], [1, 0], [0, 1]]] .* 6.67
    tri_moment = Subzero.calc_moment_inertia(tri_coords, 0.5)
    @test isapprox(tri_moment, 151743.437, atol = 0.001)
end