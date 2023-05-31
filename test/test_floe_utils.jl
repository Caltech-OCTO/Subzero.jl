@testset "Floe Utils" begin
    ext = [[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    hole1 = [[0.2, 0.3], [0.2, 0.2], [0.3, 0.2], [0.3, 0.3], [0.2, 0.3]]
    hole2 = [[0.5, 0.6], [0.5, 0.5], [0.6, 0.5], [0.6, 0.6], [0.5, 0.6]]
    poly_nohole = LG.Polygon([ext])
    poly_hole1 = LG.Polygon([ext, hole1])
    poly_hole2 = LG.Polygon([ext, hole1, hole2])
    multipoly_hole1 = LG.MultiPolygon([poly_nohole, poly_hole1])
    multipoly_hole2 = LG.MultiPolygon([poly_nohole, poly_hole1, poly_hole2])

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
    @test Subzero.hashole(multipoly_hole1)

    # Test removing holes from polygons and multipolygons
    @test [ext] == Subzero.rmholes([ext])
    @test [ext] == Subzero.rmholes([ext, hole1])
    copy_holes = [ext, hole1]
    Subzero.rmholes!(copy_holes)
    @test copy_holes == [ext]
    @test LG.equals(Subzero.rmholes(poly_nohole), poly_nohole)
    @test LG.equals(Subzero.rmholes(poly_hole1), poly_nohole)
    @test LG.equals(Subzero.rmholes(poly_hole2), poly_nohole)
    @test !Subzero.hashole(Subzero.rmholes(multipoly_hole1))

    # Test sorting regions from polygons and multipolygons
    @test Subzero.sortregions(poly_nohole) == [poly_nohole]
    poly_lst = Subzero.sortregions(multipoly_hole2)
    @test LG.equals(poly_lst[1], poly_nohole)
    @test LG.equals(poly_lst[3], poly_hole2)

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

    # Test scaling polygons
    centered_coords = [[[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], 
                        [-1.0, 1.0], [-1.0, -1.0]]]
    centered_poly = LG.Polygon(centered_coords)
    @test LG.equals(Subzero.scale(centered_poly, 1), centered_poly)
    @test LG.equals(Subzero.scale(centered_poly, 2.0),
        LG.Polygon(centered_coords .* 2.0))

    # Test seperating x and y coordiantes from PolyVec form
    x_ext, y_ext = Subzero.separate_xy([ext])
    @test x_ext == [0.0, 0.0, 1.0, 1.0, 0.0]
    @test y_ext == [1.0, 0.0, 0.0, 1.0, 1.0]

    # Test moment of intertia calculations - compared to values output my MATLAB
    poly_moment = Subzero.calc_moment_inertia([ext], [0.5, 0.5], 0.25)
    @test isapprox(poly_moment, 38.333, atol = 0.001)
    @test Subzero.calc_moment_inertia(poly_nohole, 0.25) == poly_moment
    tri_coords = [[[0, 1], [0, 0], [1, 0], [0, 1]]] .* 6.67
    tri_moment = Subzero.calc_moment_inertia(LG.Polygon(tri_coords), 0.5)
    @test isapprox(tri_moment, 50581.145, atol = 0.001)

    # Test orient_coords
    c1 = [
        [9.75e4, 7e4],
        [9.75e4, 5e4],
        [9.75e4, 5e4],
        [10.05e4, 5e4],
        [10.05e4, 7e4],
    ]
    Subzero.orient_coords!(c1)
    @test c1 == [
        [97500.0, 50000.0],
        [97500.0, 70000.0],
        [100500.0, 70000.0],
        [100500.0, 50000.0],
        [97500.0, 50000.0],
    ]
    c2 = [
        [6.5e4, 6.5e4],
        [8.5e4, 6.5e4],
        [8.5e4, 4.5e4],
        [8.5e4, 4.5e4],
        [6.5e4, 4.5e4],
    ]
    Subzero.orient_coords!(c2)
    @test c2 == [
        [6.5e4, 4.5e4],
        [6.5e4, 6.5e4],
        [8.5e4, 6.5e4],
        [8.5e4, 4.5e4],
        [6.5e4, 4.5e4],
    ]

    # Test polygon angles - some basic shapes and then compared to MATLAB
    rect_coords = [[[1.0, 1.0], [1.0, 2.0], [2.0, 2.0], [2.0, 1.0]]]
    @test Subzero.calc_poly_angles(rect_coords) == [90.0, 90.0, 90.0, 90.0]
    tri_coords = [[[0.0, 0.0], [0.0, 4.0], [3.0, 0.0]]]
    @test prod(isapprox.(
        Subzero.calc_poly_angles(tri_coords),
        [90.0, 36.8699, 53.1301],
        atol = 0.001,
    ))
    concave_tri_coords = [[[-3.0, -2.0], [0.0,0.0], [5.0, 0.0]]]
    @test prod(isapprox.(
        Subzero.calc_poly_angles(concave_tri_coords),
        [19.6538, 146.3099, 14.0362],
        atol = 0.001,
    ))
    # generate list of random polygons
    polygon_lst = voronoicells(
        rand(10),
        rand(10),
        Rectangle(Point2(0.0, 0.0), Point2(1.0, 1.0)),
    ).Cells
    for poly in polygon_lst
        @test isapprox(
            sum(Subzero.calc_poly_angles([Vector{Vector{Float64}}(poly)])),
            180 * (length(poly) - 2),
            atol = 1e-3,
        )
    end

    # Test calc_point_poly_dist - some basic shapes and compared to  MATLAB
    xpoints = [0.0, 1.1, 1.1, 2.0]
    ypoints = [0.0, 1.0, 1.1, 2.0]
    @test prod(isapprox(
        Subzero.calc_point_poly_dist(xpoints, ypoints, rect_coords),
        [sqrt(2); 0; -0.1; 0.0],
        atol = 0.001),
    )
    @test_throws AssertionError Subzero.calc_point_poly_dist(
        xpoints,
        ypoints,
        [[[1.0, 1.0], [1.0, 2.0]]],
    )
    @test Subzero.calc_point_poly_dist(
        Float64[],
        Float64[],
        rect_coords,
    ) == Float64[]
    @test_throws AssertionError Subzero.calc_point_poly_dist(
        Float64[],
        [1.0, 2.0],
        rect_coords,
    )

    # ------------------------- Test intersection of lines ---------------------
    l1 = [[[0.0, 0.0], [2.5, 0.0], [5.0, 0.0]]]
    l2 = [[[2.0, -3.0], [3.0, 0.0], [4.0, 3.0]]]
    @test Subzero.intersect_lines(l1, l2) == [3.0 0.0]
    l1 = [[[0., -1], [1, 1], [2, -1], [3, 1]]]
    l2 = [[[0., 0], [1, 0], [3, 0]]]
    @test Subzero.intersect_lines(l1, l2) == [0.5 0; 1.5 0; 2.5 0]
    l2 = [[[10., 10]]]
    @test Subzero.intersect_lines(l1, l2) == zeros(2,0)

    # -------------- Test cutting polygon through horizontal line --------------
    # Cut a hexagon through the line y = -1
    poly_coords = [[
        [2.0, -3.0],
        [0.0, 0.0],
        [2.0, 2.0],
        [6.0, 2.0],
        [8.0, 0.0],
        [6.0, -3.0],
        [2.0, -3.0,],
    ]]
    poly = LG.Polygon(Subzero.cut_polygon_coords(poly_coords, -1)[1])
    @test LG.isValid(poly)
    cut_coords = [[
        [2.0, -3.0],
        [2-4/3, -1.0],
        [22/3, -1.0],
        [6.0, -3.0],
        [2.0, -3.0],
    ]]
    @test Set(cut_coords[1]) == Set(LG.GeoInterface.coordinates(poly)[1])
    @test isapprox(
        LG.GeoInterface.coordinates(LG.centroid(poly)),
        LG.GeoInterface.coordinates(LG.centroid(LG.Polygon(cut_coords))),
        atol = 1e-6,
    )

    # Cut a c-shaped polygon through the line y = 5, creating two polygons
    poly_coords = [[
        [0.0, 0.0],
        [0.0, 10.0],
        [10.0, 10.0],
        [10.0, 0.0],
        [4.0, 0.0],
        [4.0, 6.0],
        [2.0, 6.0],
        [2.0, 0.0],
        [0.0, 0.0],
    ]]
    poly1_coords, poly2_coords = Subzero.cut_polygon_coords(poly_coords, 5)
    poly1 = LG.Polygon(poly1_coords)
    poly2 = LG.Polygon(poly2_coords)
    @test LG.isValid(poly1)
    @test LG.isValid(poly2)
    @test Set([
        [4.0, 0.0],
        [4.0, 5.0],
        [10.0, 5.0],
        [10.0, 0.0],
        [4.0, 0.0],
    ]) == Set(LG.GeoInterface.coordinates(poly1)[1])
    @test Set([
        [0.0, 0.0],
        [0.0, 5.0],
        [2.0, 5.0],
        [2.0, 0.0],
        [0.0, 0.0],
    ]) == Set(LG.GeoInterface.coordinates(poly2)[1])

    # Cut a triangle through the line y = -10 so only tip is left -> No polygon
    poly_coords = [[[0.0, 0.0], [10.0, 0.0], [5.0, -10.0], [0.0, 0.0]]]
    poly = Subzero.cut_polygon_coords(poly_coords, -10)
    @test isempty(poly)

    # Cut a rectangle through line y = 0 so only bottom edge is left -> No polygon
    poly_coords = [[
        [0.0, 0.0],
        [0.0, 5.0],
        [10.0, 5.0],
        [10.0, 0.0],
        [0.0, 0.0],
    ]]
    poly = Subzero.cut_polygon_coords(poly_coords, 0)
    @test isempty(poly)

    # Cut a rectangle through the line y = -10 so no polygon is left
    poly = Subzero.cut_polygon_coords(poly_coords, -10)
    @test isempty(poly)

    # ------------------ Test polygon splitting through hole ------------------
    # Polygon with no holes -> returns original polygon and an empty list
    poly_coords = [[
        [0.0, 0.0],
        [0.0, 5.0],
        [10.0, 5.0],
        [10.0, 0.0],
        [0.0, 0.0],
    ]]
    poly_bottom, poly_top = Subzero.split_polygon_hole(LG.Polygon(poly_coords))
    @test length(poly_bottom) == 1
    @test length(poly_top) == 0
    @test LG.area(LG.difference(poly_bottom[1], LG.Polygon(poly_coords))) == 0
    # Polygon with one hole -> creates two polygons
    poly_coords = [
        [
            [0.0, 0.0],
            [0.0, 5.0],
            [10.0, 5.0],
            [10.0, 0.0],
            [0.0, 0.0],
        ],
        [ # hole!
            [3.0, 2.0],
            [3.0, 4.0],
            [6.0, 4.0],
            [6.0, 2.0],
            [3.0, 2.0],
        ],
    ]
    poly_bottom, poly_top = Subzero.split_polygon_hole(LG.Polygon(poly_coords))
    poly_bottom_coords = [[
        [0.0, 0.0],
        [0.0, 3.0],
        [3.0, 3.0],
        [3.0, 2.0],
        [6.0, 2.0],
        [6.0, 3.0],
        [10.0, 3.0],
        [10.0, 0.0],
        [0.0, 0.0],
    ]]
    poly_top_coords = [[
        [0.0, 3.0],
        [0.0, 5.0],
        [10.0, 5.0],
        [10.0, 3.0],
        [6.0, 3.0],
        [6.0, 4.0],
        [3.0, 4.0],
        [3.0, 3.0],
        [0.0, 3.0],
    ]]
    @test length(poly_bottom) == length(poly_top) == 1
    @test LG.area(LG.difference(
        poly_bottom[1],
        LG.Polygon(poly_bottom_coords),
    )) == 0
    @test LG.area(LG.difference(poly_top[1], LG.Polygon(poly_top_coords))) == 0
    # Polygon with one holes -> creates three polygons 
    poly_coords = [
        [
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [10.0, 0.0],
            [4.0, 0.0],
            [4.0, 6.0],
            [2.0, 6.0],
            [2.0, 0.0],
            [0.0, 0.0],
        ], 
        [
            [6.0, 4.0],
            [6.0, 6.0],
            [7.0, 6.0],
            [7.0, 4.0],
            [6.0, 4.0],
        ],
    ]
    poly_bottom, poly_top = Subzero.split_polygon_hole(LG.Polygon(poly_coords))
    poly_bottom1_coords = [[
        [4.0, 0.0],
        [4.0, 5.0],
        [6.0, 5.0],
        [6.0, 4.0],
        [7.0, 4.0],
        [7.0, 5.0],
        [10.0, 5.0],
        [10.0, 0.0],
        [4.0, 0.0],
    ]]
    poly_bottom2_coords = [[
        [0.0, 0.0],
        [0.0, 5.0],
        [2.0, 5.0],
        [2.0, 0.0],
        [0.0, 0.0],
    ]]
    poly_top_coords = [[
        [0.0, 5.0],
        [0.0, 10.0],
        [10.0, 10.0],
        [10.0, 5.0],
        [7.0, 5.0],
        [7.0, 6.0],
        [6.0, 6.0],
        [6.0, 5.0],
        [4.0, 5.0],
        [4.0, 6.0],
        [2.0, 6.0],
        [2.0, 5.0],
        [0.0, 5.0],
    ]]
    @test length(poly_top) == 1
    @test length(poly_bottom) == 2
    @test LG.area(LG.difference(
        poly_bottom[1],
        LG.Polygon(poly_bottom1_coords),
    )) == 0
    @test LG.area(LG.difference(
        poly_bottom[2],
        LG.Polygon(poly_bottom2_coords),
    )) == 0
    @test LG.area(LG.difference(
        poly_top[1],
        LG.Polygon(poly_top_coords),
    )) == 0

    # Test points_in_poly for polygons
    points = [-5 5; 2.5 2.5; 2 7; 5 6; 5.5 5.5; 9 2]
    # No holes
    coords = [[[0.0, 10.0], [0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]]
    @test all(Subzero.points_in_poly(points, coords) .== 
        [false, true, true, true, true, true])
    # two holes
    coords = [
        [[0.0, 10.0], [0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]],
        [[2.0, 3], [2, 2], [3, 2], [3, 3], [2, 3]],
        [[5.0, 6], [5, 5], [6, 5], [6, 6], [5, 6]],
    ]
    @test all(Subzero.points_in_poly(points, coords) .== 
        [false, false, true, false, false, true])

    # Test points_in_poly for multi-polygons
    multi_coords = [coords, [[[12.0, 5], [15, 10], [18, 5], [12, 5]]]]
    points = vcat(points, [15 7; 11 10])
    @test all(Subzero.points_in_poly(points, multi_coords) .== 
    [false, false, true, false, false, true, true, false])

    # No coords / no points
    @test all(Subzero.points_in_poly(points, [[[Vector{Float64}()]]]) .== 
    [false, false, false, false, false, false, false, false])
    @test all(Subzero.points_in_poly([], multi_coords) .== Vector{Bool}())

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