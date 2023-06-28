@testset "Floe" begin
    FT = Float64
    # test generating monte carlo points
    file = jldopen("inputs/floe_shapes.jld2", "r")
    floe_coords = file["floe_vertices"][1:end]
    close(file)
    poly1 = LG.Polygon(Subzero.valid_polyvec!(floe_coords[1]))
    centroid1 = LG.GeoInterface.coordinates(LG.centroid(poly1))
    origin_coords = Subzero.translate(
        floe_coords[1],
        -centroid1[1],
        -centroid1[2],
    )
    xo, yo = Subzero.separate_xy(origin_coords)
    rmax = sqrt(maximum([sum(xo[i]^2 + yo[i]^2) for i in eachindex(xo)]))
    area = LG.area(poly1)
    mc_x, mc_y, status = Subzero.generate_mc_points(
        Float64,
        1000,
        origin_coords,
        rmax,
        area,
        Subzero.Status(),
        Xoshiro(1)
    )
    @test length(mc_x) == length(mc_y) && length(mc_x) > 0
    in_on = inpoly2(hcat(mc_x, mc_y), hcat(xo, yo))
    mc_in = in_on[:, 1] .|  in_on[:, 2]
    @test all(mc_in)
    @test abs(sum(mc_in)/1000 * 4 * rmax^2 - area)/area < 0.1
    @test status.tag == Subzero.active
    # Test that random number generator is working
    mc_x2, mc_y2, status2 = Subzero.generate_mc_points(
        Float64,
        1000,
        origin_coords,
        rmax,
        area,
        Subzero.Status(),
        Xoshiro(1)
    )
    @test all(mc_x .== mc_x2)
    @test all(mc_y .== mc_y2)
    @test status2.tag == Subzero.active

    mc_x3, mc_y3, status3 = Subzero.generate_mc_points(
        Float32,
        1000,
        origin_coords,
        rmax,
        area,
        Subzero.Status(),
        Xoshiro(1)
    )
    @test status3.tag == Subzero.active
    @test eltype(mc_x3) == eltype(mc_y3) == Float32

    rhombus_coords = [[
        [0.0, 0.0],
        [40.0, 50.0],
        [100.0, 50.0],
        [80.0, 0.0],
        [0.0, 0.0],
    ]]
    rhombus_centroid = [54.285, 23.809]

    small_grid = RegRectilinearGrid(
        (-150, 150),
        (-150, 150),
        5,   # Δx
        10,  # Δy
    )

    medium_grid = RegRectilinearGrid(
        (-150, 150),
        (-150, 150),
        40,   # Δx
        40,  # Δy
    )

    large_grid = RegRectilinearGrid(
        (-150, 150),
        (-150, 150),
        100,   # Δx
        100,  # Δy
    )

    x, y = Subzero.generate_floe_points(
        rhombus_coords,
        rhombus_centroid,
        1,  # number of points per grid cell
        small_grid,
    )

    # Test InteractionFields enum
    interactions = range(1, 7)'
    @test interactions[floeidx] == 1
    @test interactions[xpoint] == 4
    @test interactions[overlap] == 7

    # Test with coords inputs
    floe_from_coords = Floe(
        floe_coords[1],
        0.5,
        0.01,
        u = 0.2,
        rng = Xoshiro(1),
    )
    @test typeof(floe_from_coords) <: Floe
    @test floe_from_coords.u == 0.2
    @test 0.49 <= floe_from_coords.height <= 0.51
    @test floe_from_coords.centroid == centroid1
    @test floe_from_coords.area == area
    @test floe_from_coords.status.tag == Subzero.active
    
    # Test with polygon input
    floe_from_poly = Floe(
        poly1,
        0.5,
        0.01,
        v = -0.2,
        rng = Xoshiro(1),
    )
    @test typeof(floe_from_poly) <: Floe
    @test floe_from_poly.u == 0.0
    @test floe_from_poly.v == -0.2
    @test 0.49 <= floe_from_poly.height <= 0.51
    @test floe_from_poly.centroid == centroid1
    @test floe_from_poly.area == area
    @test floe_from_poly.status.tag == Subzero.active

    # Test poly_to_floes
    rect_poly = [[[0.0, 0.0], [0.0, 5.0], [10.0, 5.0], [10.0, 0.0], [0.0, 0.0]]]

    c_poly_hole = [[[0.0, 0.0], [0.0, 10.0], [10.0, 10.0], [10.0, 0.0],
        [4.0, 0.0], [4.0, 6.0], [2.0, 6.0], [2.0, 0.0], [0.0, 0.0]],
        [[6.0, 4.0], [6.0, 6.0], [7.0, 6.0], [7.0, 4.0], [6.0, 4.0]]]

    # Test polygon with no holes
    floe_arr = Subzero.poly_to_floes(
        FT,
        LG.Polygon(rect_poly),
        0.25,
        0.01,
    )
    @test length(floe_arr) == 1
    @test typeof(floe_arr[1]) <: Floe
    @test !Subzero.hashole(floe_arr.coords[1])

    # Test with polygon below minimum floe area
    floe_arr = Subzero.poly_to_floes(
        FT,
        LG.Polygon(rect_poly),
        0.25,
        0.01,
        min_floe_area = 55
    )
    @test isempty(floe_arr)

    # Test with polygon with a hole that is split into 3 polyons
    floe_arr = Subzero.poly_to_floes(
        FT,
        LG.Polygon(c_poly_hole),
        0.25,
        0.01,
    )
    @test length(floe_arr) == 3
    @test !any(Subzero.hashole.(floe_arr.coords))
    @test typeof(floe_arr) <: StructArray{<:Floe}

    # Test with multipolygon 
    floe_arr = Subzero.poly_to_floes(
        FT,
        LG.MultiPolygon([c_poly_hole, rect_poly]),
        0.25,
        0.01
    )
    @test length(floe_arr) == 4
    @test typeof(floe_arr) <: StructArray{<:Floe}

    # Test with multipolygon with minimum area
    floe_arr = Subzero.poly_to_floes(
        FT,
        LG.MultiPolygon([c_poly_hole, rect_poly]),
        0.25,
        0.01,
        min_floe_area = 30
    )
    @test length(floe_arr) == 2
    @test typeof(floe_arr) <: StructArray{<:Floe}

    # Test initialize_floe_field from coord list
    grid = RegRectilinearGrid(
        (-8e4, 8e4),
        (-8e4, 8e4),
        1e4,
        1e4,
    )
    nbound = CollisionBoundary(North, grid)
    sbound = CollisionBoundary(South, grid)
    ebound = CollisionBoundary(East, grid)
    wbound = CollisionBoundary(West, grid)
    domain_no_topo = Domain(nbound, sbound, ebound, wbound)
    island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
    topo1 = [[[-8e4, -8e4], [-8e4, 8e4], [-6e4, 8e4], [-5e4, 4e4], [-6e4, -8e4], [-8e4, -8e4]]]
    domain_with_topo = Domain(
        nbound,
        sbound,
        ebound,
        wbound,
        StructVector([TopographyElement(t) for t in [island, topo1]])
    )
    topo_polys = LG.MultiPolygon([island, topo1])

    # From file without topography -> floes below recommended area
    floe_arr = (@test_logs (:warn, "Some user input floe areas are less than the suggested minimum floe area.") initialize_floe_field(
        FT,
        floe_coords,
        domain_no_topo,
        0.5,
        0.1,
    ))
    nfloes = length(floe_coords)
    @test typeof(floe_arr) <: StructArray{<:Floe}
    @test length(floe_arr) == nfloes
    @test all(floe_arr.id .== range(1, nfloes))

    # From file with small domain -> floes outside of domain
    small_grid = RegRectilinearGrid(
        (-5e4, 5e4),
        (-5e4, 5e4),
        1e4,
        1e4,
    )
    small_domain_no_topo = Domain(
        CollisionBoundary(North, small_grid),
        CollisionBoundary(South, small_grid),
        CollisionBoundary(East, small_grid),
        CollisionBoundary(West, small_grid),
    )
    floe_arr = (@test_logs (:warn, "Some floe centroids are out of the domain.") initialize_floe_field(
        FT,
        floe_coords,
        small_domain_no_topo,
        0.5,
        0.1,
        min_floe_area = 1e5,
    ))
    @test typeof(floe_arr) <: StructArray{<:Floe}
    @test all(floe_arr.area .> 1e5)

    # From file with topography -> no intersection between topo and floes
    floe_arr = initialize_floe_field(
        FT,
        floe_coords,
        domain_with_topo,
        0.5,
        0.1,
        min_floe_area = 10,
        rng = Xoshiro(0)
    )
    @test typeof(floe_arr) <: StructArray{<:Floe}
    @test all([LG.area(
        LG.intersection(LG.Polygon(c), topo_polys)
    ) for c in floe_arr.coords] .< 1e-6)

    # Test generate_voronoi_coords - general case
    domain_coords = [[[[1, 2], [1.5, 3.5], [1, 5], [2.5, 5], [2.5, 2], [1, 2]]]]
    bounding_box = [[[1, 2], [1, 5], [2.5, 5], [2.5, 2], [1, 2]]]
    voronoi_coords = Subzero.generate_voronoi_coords(
        10,
        [1.5, 3],
        [1, 2],
        domain_coords,
        Xoshiro(1),
        10,
        max_tries = 20, # 20 tries makes it very likely to reach 10 polygons
    )
    bounding_poly = LG.Polygon(bounding_box)
    @test length(voronoi_coords) == 10
    for c in voronoi_coords
        fpoly = LG.Polygon(c)
        @test isapprox(
            LG.area(LG.intersection(fpoly, bounding_poly)),
            LG.area(fpoly),
            atol = 1e-3,
        )
        @test LG.isValid(fpoly)
    end

    # Test warning and no points generated
    warning_str = "Only 0 floes were able to be generated in 10 tries during \
        voronoi tesselation."
    @test @test_logs (:warn, warning_str) Subzero.generate_voronoi_coords(
        0, # Don't generate any points
        [1.5, 3],
        [1, 2],
        domain_coords,
        Xoshiro(1),
        10, # Higher min to warn so we can test warning
    ) == Vector{Vector{Vector{Float64}}}()
    
    # Test initialize_floe_field with voronoi
    floe_arr = initialize_floe_field(
        FT,
        25,
        [0.5],
        domain_with_topo,
        0.5,
        0.1,
        min_floe_area = 1e4,
        rng = Xoshiro(1)
    )
    @test isapprox(
        sum(floe_arr.area)/(1.6e5*1.6e5 - LG.area(topo_polys)),
        0.5,
        atol = 1e-1
    )
    @test all(floe_arr.area .> 1e4)
    @test all([LG.area(
        LG.intersection(LG.Polygon(c), topo_polys)
    ) for c in floe_arr.coords] .< 1e-6)
    nfloes = length(floe_arr)
    @test all(floe_arr.id .== range(1, nfloes))

    concentrations = [1 0.3; 0 0.5]
    rng = Xoshiro(2)
    floe_arr = initialize_floe_field(
        25,
        concentrations,
        domain_with_topo,
        0.5,
        0.1,
        min_floe_area = 1e4,
        rng = rng
    )
    nfloes = length(floe_arr)
    floe_polys = [LG.Polygon(f) for f in floe_arr.coords]
    first_cell = [[[-8e4, -8e4], [-8e4, 0], [0, 0], [0, -8e4], [-8e4, -8e4]]]
    for j in 1:2
        for i in 1:2
            cell = LG.Polygon(
                Subzero.translate(first_cell, 8e4*(j-1), 8e4*(i-1))
            )
            open_cell_area = LG.area(LG.difference(cell, topo_polys))
            c = concentrations[i, j]
            floes_in_cell = [LG.intersection(p, cell) for p in floe_polys]
            @test c - 100eps() <= sum(LG.area.(floes_in_cell))/open_cell_area < 1 + eps()
        end
    end
    @test all([LG.area(
        LG.intersection(p, topo_polys)
    ) for p in floe_polys] .< 1e-3)
    @test all([LG.isValid(p) for p in floe_polys])
    @test all(floe_arr.id .== range(1, nfloes))

    @test typeof(initialize_floe_field(
        Float32,
        25,
        concentrations,
        domain_with_topo,
        0.5,
        0.1,
    )) <: StructArray{<:Floe{Float32}}
end