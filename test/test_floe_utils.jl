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
    c1_new = Subzero.orient_coords(c1)
    @test c1_new == [
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
    c2_new = Subzero.orient_coords(c2)
    @test c2_new == [
        [6.5e4, 4.5e4],
        [6.5e4, 6.5e4],
        [8.5e4, 6.5e4],
        [8.5e4, 4.5e4],
        [6.5e4, 4.5e4],
    ]

    # Test calc_point_poly_dist - some basic shapes and compared to  MATLAB
    rect_coords = [[[1.0, 1.0], [1.0, 2.0], [2.0, 2.0], [2.0, 1.0]]]
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
    @test issetequal(
        Subzero.intersect_lines(l1, l2),
        Set([(3.0, 0.0)]),
    )
    l1 = [[[0., -1], [1, 1], [2, -1], [3, 1]]]
    l2 = [[[0., 0], [1, 0], [3, 0]]]
    @test issetequal(
        Subzero.intersect_lines(l1, l2),
        Set([(0.5, -0.0), (1.5, 0), (2.5, -0.0)]),
    )
    l2 = [[[10., 10]]]
    @test issetequal(
        Subzero.intersect_lines(l1, l2),
        Set{Tuple{Float64, Float64}}(),
    )

    # ------------------------- Test finding shared points ---------------------
    two_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[20.0, 0.0], [20.0, 20.0], [40.0, 20.0], [40.0, 0.0], [20.0, 0.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        two_shared_v[1][1],
        two_shared_v[2],
    ) == [1, 2]

    @test Subzero.which_points_on_edges(
        two_shared_v[1][1],
        two_shared_v[2],
    ) == [3, 4]

    three_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 10.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[40.0, 20.0], [40.0, 0.0], [20.0, 0.0], [20.0, 10.0], [20.0, 20.0], [40.0, 20.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        three_shared_v[1][1],
        three_shared_v[2],
    ) == [3, 4, 5]

    @test Subzero.which_points_on_edges(
        three_shared_v[1][1],
        three_shared_v[2],
    ) == [3, 4, 5]

    deleteat!(three_shared_v[2][1], 4)
    @test Subzero.which_points_on_edges(
        three_shared_v[1][1],
        three_shared_v[2],
    ) == [3, 4, 5]

    four_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 18.0], [20.0, 15.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[20.0, 18.0], [20.0, 20.0], [40.0, 20.0], [40.0, 0.0], [20.0, 0.0], [20.0, 15.0], [20.0, 18.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        four_shared_v[1][1],
        four_shared_v[2],
    ) == [1, 2, 5, 6]

    offset_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [20.0, 0.0], [0.0, 0.0]]],
        [[[5.0, 20.0], [5.0, 25.0], [25.0, 25.0], [25.0, 20.0], [10.0, 20.0], [5.0, 20.0]]]
    ]
    @test Subzero.which_points_on_edges(
        offset_shared_v[1][1],
        offset_shared_v[2],
    ) == [3]
    @test Subzero.which_points_on_edges(
        offset_shared_v[2][1],
        offset_shared_v[1],
    ) == [1, 5]

    triange_shared_v = [
        [[[0.0, 0.0], [0.0, 20.0], [20.0, 20.0], [5.0, 5.0], [0.0, 0.0]]],
        [[[0.0, 0.0], [5.0, 5.0], [20.0, 20.0], [20.0, 0.0], [0.0, 0.0]]]
    ]
    @test Subzero.which_vertices_match_points(
        triange_shared_v[1][1],
        triange_shared_v[2],
    ) == [1, 2, 3]

    @test Subzero.which_points_on_edges(
        triange_shared_v[1][1],
        triange_shared_v[2],
    ) == [1, 3, 4]

    deleteat!(triange_shared_v[2][1], 2)
    @test Subzero.which_points_on_edges(
        triange_shared_v[1][1],
        triange_shared_v[2],
    ) == [1, 3, 4]

    squares_midpoint = (20.0, 10.0)
    triangle_midpoint = (10.0, 10.0)
    @test Subzero.find_shared_edges_midpoint(
        two_shared_v[1],
        two_shared_v[2],
    ) == squares_midpoint
    @test Subzero.find_shared_edges_midpoint(
        three_shared_v[1],
        three_shared_v[2],
    ) == squares_midpoint
    @test Subzero.find_shared_edges_midpoint(
        four_shared_v[1],
        four_shared_v[2],
    ) == squares_midpoint
    @test Subzero.find_shared_edges_midpoint(
        triange_shared_v[1],
        triange_shared_v[2],
    ) == (10.0, 10.0)
    @test Subzero.find_shared_edges_midpoint(
        offset_shared_v[2],
        offset_shared_v[1],
    ) == (7.5, 20)

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