using Test, Subzero

    # Boundaries using BoundaryCoords
    FT = Float64
    g = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
    b1 = Subzero.PeriodicBoundary(North, g)
    b2 = Subzero.OpenBoundary(East, g)
    b3 = Subzero.CollisionBoundary(West, g)
    b4 = Subzero.PeriodicBoundary(South, g)
    b5 = Subzero.MovingBoundary(South, g, 1.0, 2.0)
    @test b1.val == 3e5
    @test typeof(b1) == Subzero.PeriodicBoundary{North, Float64}
    @test b1.coords == [[[-2e5, 3e5], [-2e5, 4.5e5], [6e5, 4.5e5], [6e5, 3e5], [-2e5, 3e5]]]
    @test b2.val == 4e5
    @test typeof(b2) == Subzero.OpenBoundary{East, Float64}
    @test b2.coords == [[[4e5, -1.5e5], [4e5, 4.5e5], [6e5, 4.5e5], [6e5, -1.5e5], [4e5, -1.5e5]]]
    @test b3.val == 0.0
    @test typeof(b3) == Subzero.CollisionBoundary{West, Float64}
    @test b3.coords == [[[-2e5, -1.5e5], [-2e5, 4.5e5], [0.0, 4.5e5], [0.0, -1.5e5], [-2e5, -1.5e5]]]
    @test b4.val == 0.0
    @test typeof(b4) == Subzero.PeriodicBoundary{South, Float64}
    @test b4.coords == [[[-2e5, -1.5e5], [-2e5, 0.0], [6e5, 0.0], [6e5, -1.5e5], [-2e5, -1.5e5]]]
    @test b5.u == 1.0
    @test b5.v == 2.0

    # Test changing CollisionBoundary values and can't change other boundary types
    b5.val = 1.0
    @test b5.val == 1.0
    @test_throws Exception b4.val = 1.0

    # Creation of Float32 and Float64 Boundary
    b32 = Subzero.OpenBoundary(Float32, North, g)
    @test typeof(b32.val) == Float32
    @test typeof(b32) == Subzero.OpenBoundary{North, Float32}
    b64 = Subzero.OpenBoundary(Float64, North, g)
    @test typeof(b64.val) == Float64
    @test typeof(b64) == Subzero.OpenBoundary{North, Float64}

    # Periodic Compat
    @test !Subzero.periodic_compat(b1, b2)
    @test !Subzero.periodic_compat(b2, b1)
    @test Subzero.periodic_compat(b1, b4)
    @test Subzero.periodic_compat(b2, b3)


    # Test compression boundaries movement - TODO: Move these to different file
    nc_boundary = MovingBoundary(North, grid, 0.0, -0.1)
    nc_coords = deepcopy(nc_boundary.coords)
    sc_boundary = MovingBoundary(South, grid, 0.0, 0.1)
    sc_coords = deepcopy(sc_boundary.coords)
    ec_boundary = MovingBoundary(East, grid, 0.1, 0.0)
    ec_coords = deepcopy(ec_boundary.coords)
    wc_boundary = MovingBoundary(West, grid, 0.1, 0.0)
    wc_coords = deepcopy(wc_boundary.coords)
    cdomain = Domain(nc_boundary, sc_boundary, ec_boundary, wc_boundary)
    Subzero.update_boundaries!(cdomain, 10)
    Subzero.translate!(nc_coords, 0, -1)
    @test nc_coords == nc_boundary.coords
    @test nc_boundary.val == 1e5 - 1
    Subzero.translate!(sc_coords, 0, 1)
    @test sc_coords == sc_boundary.coords
    @test sc_boundary.val == -1e5 + 1
    Subzero.translate!(ec_coords, 1, 0)
    @test ec_coords == ec_boundary.coords
    @test ec_boundary.val == 1e5 + 1
    Subzero.translate!(wc_coords, 1, 0)
    @test wc_coords == wc_boundary.coords
    @test wc_boundary.val == -1e5 + 1