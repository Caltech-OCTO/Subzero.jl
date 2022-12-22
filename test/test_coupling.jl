@testset "Coupling" begin
    grid = Subzero.RegRectilinearGrid(-10, 10, -8, 8, 2, 4)
    # Test find_cell_indices
    @test Subzero.find_cell_indices([], [], grid) == (Vector{Int}(undef, 0), Vector{Int}(undef, 0))
    @test Subzero.find_cell_indices([-10.5, -10, -10, -6.5, -6, -4, 10, 10.5, 12],
                                    [0.0, 6.0, -8.0, 4.5, 0.0, 5.0, -8.0, 0.0, 0.0], grid) ==
        ([1, 1, 1, 3, 3, 4, 11, 11, 12], [3, 5, 1, 4, 3, 4, 1, 3, 3])

    # Test filter_oob_points
    open_bound = Subzero.OpenBoundary(grid, Subzero.East())
    periodic_bound = Subzero.PeriodicBoundary(grid, Subzero.East())
    p = [-12 -10 -8 -6 0 4 4 10 12 12; 5 -6 4 10 -10 8 -8 -6 4 10]
    x = p[1, :]
    y = p[2, :]
    # all bounds non-periodic - points outside of grid in x and y
    @test Subzero.filter_oob_points(p, x, y, grid, open_bound, open_bound) == 
        (p[:, (-10 .<= x .<= 10) .& (-8 .<= y .<= 8)],
        x[(-10 .<= x .<= 10) .& (-8 .<= y .<= 8)],
        y[(-10 .<= x .<= 10) .& (-8 .<= y .<= 8)])
    warning_str = "A floe longer than the domain passed through a periodic boundary. It was removed to prevent overlap."
    # y bounds periodic - points outside of grid in south and north so all points removed
    @test (@test_logs (:warn, warning_str) Subzero.filter_oob_points(p, x, y, grid, open_bound, periodic_bound)) ==
        (Matrix{Int}(undef, 2, 0), Vector{Int}(undef, 0), Vector{Int}(undef, 0))
    # x bounds periodic - points outside of grid in east and west so all points removed
    @test (@test_logs (:warn, warning_str) Subzero.filter_oob_points(p, x, y, grid, periodic_bound, open_bound)) == 
        (Matrix{Int}(undef, 2, 0), Vector{Int}(undef, 0), Vector{Int}(undef, 0))
    # all bounds periodic - points outside of grid in all 4 directions so all points removed
    @test (@test_logs (:warn, warning_str) Subzero.filter_oob_points(p, x, y, grid, periodic_bound, periodic_bound)) == 
        (Matrix{Int}(undef, 2, 0), Vector{Int}(undef, 0), Vector{Int}(undef, 0))
    # y bounds periodic with points only outside of periodic in north direction - filter open bounds points
    pn =  [-12 -10 -8 -6 4 10 12 12; 5 -6 4 10 8 -6 4 10]
    xn = pn[1, :]
    yn = pn[2, :]
    @test Subzero.filter_oob_points(pn, pn[1, :], pn[2, :], grid, open_bound, periodic_bound) ==
        (pn[:, -10 .<= xn .<= 10],
         xn[-10 .<= xn .<= 10],
         yn[-10 .<= xn .<= 10])
    # x bounds periodic with points only outside of periodic in east direction - filter open bounds points
    pe = [-8 -6 0 4 4 10 12 12; 4 10 -10 8 -8 -6 4 10]
    xe = pe[1, :]
    ye = pe[2, :]
    @test Subzero.filter_oob_points(pe[:,1:end .!= 2], xe[1:end .!= 2], ye[1:end .!= 2], grid, periodic_bound, open_bound) ==
        (pe[:, -8 .<= ye .<= 8],
         xe[-8 .<= ye .<= 8],
         ye[-8 .<= ye .<= 8])
    # all bounds periodic - points outside of grid in only north and east directiom so no points removed
    pne = [-8 -6 4  10 12 12; 4 10 8 -6 4 10]
    xne = pne[1, :]
    yne = pne[2, :]
    @test Subzero.filter_oob_points(pne, xne, yne, grid, periodic_bound, periodic_bound) == (pne, xne, yne)

    # Test find_interp_knots
    xg = 0:10:80
    # in bounds
    @test Subzero.find_interp_knots([4], xg, 2, open_bound) == (0:10:60, 1:7)
    @test Subzero.find_interp_knots([4], xg, 2, periodic_bound) == (0:10:60, 1:7)
    # out of bounds to left
    @test Subzero.find_interp_knots([0, 1], xg, 2, open_bound) == (0:10:30, 1:4)
    @test Subzero.find_interp_knots([0, 1], xg, 2, periodic_bound) == (-40:10:30, [5:8; 1:4])
    # out of bounds to right
    @test Subzero.find_interp_knots([8, 9], xg, 1, open_bound) == (50:10:80, 6:9)
    @test Subzero.find_interp_knots([8, 9], xg, 1, periodic_bound) == (50:10:100, [6:8; 1:3])
    # buffer out of bounds on both sides
    @test Subzero.find_interp_knots(1:8, xg, 2, open_bound) == (0:10:80, 1:9)
    @test Subzero.find_interp_knots(1:8, xg, 2, periodic_bound) == (-30:10:100, [6:8; 1:8; 1:3])
    # floe out of bounds on both sides - would be filtered out for periodic in advance
    @test Subzero.find_interp_knots(0:9, xg, 2, open_bound) == (0:10:80, 1:9)

    #Test cell_coords
    cell = Subzero.center_cell_coords(2, 3, grid, periodic_bound, periodic_bound)
    cell_poly = LibGEOS.Polygon(cell)
    @test LibGEOS.area(cell_poly)::Float64 == 8
    @test LibGEOS.GeoInterface.coordinates(cell_poly) == 
        [[[-9, -2], [-9, 2], [-7, 2], [-7, -2], [-9, -2]]]
    @test Subzero.center_cell_coords(1, 1, grid, open_bound, open_bound) ==
        [[[-10, -8], [-10, -6], [-9, -6], [-9, -8], [-10, -8]]]
    @test Subzero.center_cell_coords(11, 6, grid,  periodic_bound, periodic_bound) ==
        [[[9, 10], [9, 14], [11, 14], [11, 10], [9, 10]]]
    @test Subzero.center_cell_coords(11, 6, grid,  open_bound, open_bound) == 
        [[[9, 8], [9, 8], [10, 8], [10, 8], [9, 8]]]
    @test Subzero.center_cell_coords(11, 6, grid,  open_bound, periodic_bound) ==
        [[[9, 8], [9, 8], [11, 8], [11, 8], [9, 8]]]
    @test Subzero.center_cell_coords(11, 6, grid,  periodic_bound, open_bound) ==
        [[[9, 10], [9, 14], [10, 14], [10, 10], [9, 10]]]
    

    # Test aggragate_grid_stress!
    ocean = Subzero.Ocean(grid, 0, 0, 0)
    floe1 = Subzero.Floe([[[1.0, 2], [1, 6], [3, 6], [3, 2], [1, 2]]], 0.5, 0.0)
    # floe is only in cell (7,4) so others will not contribute due to lack of area
    Subzero.aggragate_grid_force!([6, 7, 6, 7, 6, 7], [4, 4, 3, 3, 4, 4], ones(6), 2ones(6), floe1,
                                  ocean, grid, open_bound, open_bound)
    @test ocean.fx[4, 6] == ocean.fx[3, 6] == ocean.fx[3, 7] == 0
    @test ocean.fx[4,7] == 8
    @test ocean.fy[4, 6] == ocean.fy[3, 6] == ocean.fy[3, 7] == 0
    @test ocean.fy[4, 7] == 16
    @test ocean.si_area[4, 6] == ocean.si_area[3, 6] == ocean.si_area[3, 7] == 0
    @test ocean.si_area[4, 7] == 8 

    floe2 = Subzero.Floe([[[2.0, -4], [2, 0], [6, 0], [6, -4], [2, -4]]], 0.5, 0.0)
    Subzero.aggragate_grid_force!([7, 7, 8, 8, 9, 9], [2, 3, 3, 3, 2, 2], ones(6), 2ones(6), floe2,
                                  ocean, grid, periodic_bound, periodic_bound)
    @test ocean.fx[2, 7] == ocean.fx[3, 7] == ocean.fx[2, 9]  == 2
    @test ocean.fx[3, 9] == ocean.fx[2, 8] == 0
    @test ocean.fx[3, 8] == 4
    @test ocean.fy[2, 7] == ocean.fy[3, 7] == ocean.fy[2, 9] == 4
    @test ocean.fy[3, 8] == 8
    @test ocean.si_area[2, 7] == ocean.si_area[3, 7] == ocean.si_area[2, 9] == 2
    @test ocean.si_area[3, 8] == 4
    @test ocean.si_area[3, 9] == ocean.si_area[2, 8] == 0

    ocean = Subzero.Ocean(grid, 0, 0, 0)
    floe3 = Subzero.Floe([[[7.0, 6], [7, 12], [12, 12], [12, 6], [7, 6]]], 0.5, 0.0)
    Subzero.aggragate_grid_force!([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6], ones(6), 2ones(6), floe3,
                                  ocean, grid, open_bound, open_bound)
    @test ocean.fx[5, 10] == 4
    @test ocean.fx[5, 11] == 2
    @test sum(ocean.fx) == 6

    ocean = Subzero.Ocean(grid, 0, 0, 0)
    Subzero.aggragate_grid_force!([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6], ones(6), 2ones(6), floe3,
                                  ocean, grid, periodic_bound, open_bound)

    @test ocean.fx[1, 10]  == 8
    @test ocean.fx[2, 11]  == 2
    @test ocean.fx[2, 10] == ocean.fx[1, 11] == 4
    @test sum(ocean.fx) == 18

    ocean = Subzero.Ocean(grid, 0, 0, 0)
    Subzero.aggragate_grid_force!([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6], ones(6), 2ones(6), floe3,
                                  ocean, grid, open_bound, periodic_bound)
    @test ocean.fx[5, 10] == ocean.fx[5, 1] == 4
    @test ocean.fx[5, 2] == 2
    @test sum(ocean.fx) == 10

    ocean = Subzero.Ocean(grid, 0, 0, 0)
    Subzero.aggragate_grid_force!([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6], ones(6), 2ones(6), floe3,
                                  ocean, grid, periodic_bound, periodic_bound)
    @test ocean.fx[1, 10] == ocean.fx[1, 1] == 8
    @test ocean.fx[2, 10] == ocean.fx[2, 1] == ocean.fx[1, 2] == 4
    @test ocean.fx[2, 2] == 2
    @test sum(ocean.fx) == 30

    # Test floe_OA_forcings!
end