@testset "Coupling" begin
    grid = Subzero.RegRectilinearGrid(-10, 10, -8, 8, 2, 4)
    # Test find_cell_indices
    @test Subzero.find_cell_indices([], [], grid) == (Vector{Int}(undef, 0), Vector{Int}(undef, 0))
    @test Subzero.find_cell_indices([-10, -10, -7.5, -6, -4, 10], [6, -8, 4.5, 0, 5, -8], grid) ==
        ([1, 1, 2, 3, 4, 11], [4, 1, 4, 3, 4, 1])

    # Test filter_oob_points
    open_bound = Subzero.OpenBoundary(grid, Subzero.East())
    periodic_bound = Subzero.PeriodicBoundary(grid, Subzero.East())
    p = [-12 -10 -8 -6 0 4 4 10 12 12; 5 -6 4 10 -10 8 -8 -6 4 10]
    xr = [-12, -10, -8, -6, 0, 4, 4, 10, 12, 12]
    yr = [5, -6, 4, 10, -10, 8, -8, -6, 4, 10]
    # all bounds non-periodic
    @test Subzero.filter_oob_points(p, xr, yr, grid, open_bound, open_bound) == 
        (p[:, (-10 .<= xr .<= 10) .& (-8 .<= yr .<= 8)],
        xr[(-10 .<= xr .<= 10) .& (-8 .<= yr .<= 8)],
        yr[(-10 .<= xr .<= 10) .& (-8 .<= yr .<= 8)])
    @test Subzero.filter_oob_points(p, xr, yr, grid, open_bound, periodic_bound) ==
        (p[:, -10 .<= xr .<= 10],
        xr[-10 .<= xr .<= 10],
        yr[-10 .<= xr .<= 10])
    @test Subzero.filter_oob_points(p, xr, yr, grid, periodic_bound, open_bound) == 
        (p[:, -8 .<= yr .<= 8],
        xr[-8 .<= yr .<= 8],
        yr[-8 .<= yr .<= 8])
    @test Subzero.filter_oob_points(p, xr, yr, grid, periodic_bound, periodic_bound) == 
        (p, xr, yr)

    # Test find_interp_knots
    xg = 0:10:80
    # in bounds
    @test Subzero.find_interp_knots([3, 4], xg, 2, open_bound) == (0:10:60, 1:7)
    @test Subzero.find_interp_knots([3, 4], xg, 2, periodic_bound) == (0:10:60, 1:7)
    # out of bounds to left
    @test Subzero.find_interp_knots([0, 1], xg, 2, open_bound) == (0:10:30, 1:4)
    @test Subzero.find_interp_knots([0, 1], xg, 2, periodic_bound) == (-30:10:30, [6:8; 1:4])
    # out of bounds to right
    @test Subzero.find_interp_knots([8, 9], xg, 1, open_bound) == (60:10:80, 7:9)
    @test Subzero.find_interp_knots([8, 9], xg, 1, periodic_bound) == (60:10:100, [7:9; 2:3])
    # buffer out of bounds on both sides
    @test Subzero.find_interp_knots(1:8, xg, 2, open_bound) == (0:10:80, 1:9)
    @test Subzero.find_interp_knots(1:8, xg, 2, periodic_bound) == (-20:10:100, [7:8; 1:9; 2:3])
    # floe out of bounds on both sides
    @test Subzero.find_interp_knots(0:9, xg, 2, open_bound) == (0:10:80, 1:9)
    warning_str = "A floe longer than the domain passed through a periodic boundary. It was removed to prevent overlap."
    @test (@test_logs (:warn, warning_str) Subzero.find_interp_knots(0:9, xg, 2, periodic_bound)) ==
        (Vector{Float64}(undef, 0), Vector{Float64}(undef, 0))

    #Test cell_coords
    cell = Subzero.cell_coords(2, 3, grid)
    cell_poly = LibGEOS.Polygon(cell)
    @test LibGEOS.area(cell_poly)::Float64 == 8
    @test LibGEOS.GeoInterface.coordinates(cell_poly) == 
        [[[-8, 0], [-8, 4], [-6, 4], [-6, 0], [-8, 0]]]

    # Test aggragate_grid_stress!
    ocean = Subzero.Ocean(grid, 0, 0, 0)
    floe1 = Subzero.Floe(LibGEOS.Polygon([[[1.0, 2], [1, 6], [3, 6], [3, 2], [1, 2]]]), 0.5, 0.0)
    Subzero.aggragate_grid_stress!([6, 7, 6, 7, 6, 7], [4, 4, 3, 3, 4, 4], ones(6), 2ones(6), floe1, ocean, grid)
    @test ocean.fx[4, 6] == ocean.fx[4, 7] == ocean.fx[3, 6] == ocean.fx[3, 7] == 2
    @test ocean.fy[4, 6] == ocean.fy[4, 7] == ocean.fy[3, 6] == ocean.fy[3, 7] == 4
    @test ocean.si_area[4, 6] == ocean.si_area[4, 7] == ocean.si_area[3, 6] == ocean.si_area[3, 7] == 2
end