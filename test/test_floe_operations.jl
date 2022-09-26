@testset "Floe Operations" begin
    ext = [[0.0, 1.0],  [0.0, 0.0],  [1.0, 0.0],  [1.0, 1.0], [0.0, 1.0]]
    hole1 = [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3], [0.2, 0.3]]
    hole2 = [[0.5, 0.6], [0.5, 0.5], [0.6, 0.5], [0.6, 0.6], [0.5, 0.6]]
    poly_nohole = LibGEOS.Polygon([ext])
    poly_hole1 = LibGEOS.Polygon([ext, hole1])
    poly_hole2 = LibGEOS.Polygon([ext, hole1, hole2])
    multipoly_hole1 = LibGEOS.MultiPolygon([poly_nohole, poly_hole1])
    multipoly_hole2 = LibGEOS.MultiPolygon([poly_nohole, poly_hole1,
                                                         poly_hole2])

    # Test predicate hashole for polygons and multipolygons
    @test !Subzero.hashole(poly_nohole)
    @test Subzero.hashole(poly_hole1)
    @test Subzero.hashole(poly_hole2)
    @test Subzero.hashole(multipoly_hole1)
    # Test removing holes from polygons and multipolygons
    @test LibGEOS.equals(Subzero.rmholes(poly_nohole), poly_nohole)
    @test LibGEOS.equals(Subzero.rmholes(poly_hole1), poly_nohole)
    @test LibGEOS.equals(Subzero.rmholes(poly_hole2), poly_nohole)
    @test !Subzero.hashole(Subzero.rmholes(multipoly_hole1))
    #Test sorting regions from polygons and multipolygons
    @test Subzero.sortregions(poly_nohole) == [poly_nohole]
    poly_lst = Subzero.sortregions(multipoly_hole2)
    @test LibGEOS.equals(poly_lst[1], poly_nohole)
    @test LibGEOS.equals(poly_lst[3], poly_hole2)
    #Test translating coordinates and polygons
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
    #Test scaling polygons
    centered_coords = [[[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], 
                        [-1.0, 1.0], [-1.0, -1.0]]]
    centered_poly = LibGEOS.Polygon(centered_coords)
    @test LibGEOS.equals(Subzero.scale(centered_poly, 1), centered_poly)
    @test LibGEOS.equals(Subzero.scale(centered_poly, 2.0),
        LibGEOS.Polygon(centered_coords .* 2.0))
end