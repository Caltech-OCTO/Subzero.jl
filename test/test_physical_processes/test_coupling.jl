@testset "Coupling" begin
    @testset "Coupling Helper Functions" begin
        grid = Subzero.RegRectilinearGrid(
            Float64,
            (-10, 10),
            (-8, 8),
            2,
            4,
        )
        # Test find_cell_indices
        xpoints = [-10.5, -10, -10, -6.5, -6, -4, 10, 10.5, 12]
        ypoints = [0.0, 6.0, -8.0, 4.5, 0.0, 5.0, -8.0, 0.0, 0.0]
        xidx = [1, 1, 1, 3, 3, 4, 11, 11, 12]
        yidx = [3, 5, 1, 4, 3, 4, 1, 3, 3]
        for i in eachindex(xpoints)
            @test  Subzero.find_cell_index(
                xpoints[i],
                ypoints[i],
                grid,
            ) == (xidx[i], yidx[i])
        end

        # Test filter_oob_points
        open_bound = Subzero.OpenBoundary(grid, Subzero.East())
        periodic_bound = Subzero.PeriodicBoundary(grid, Subzero.East())
        x = [-12, -10, -8, -6, 0, 4, 4, 10, 12, 12]
        y = [5, -6, 4, 10, -10, 8, -8, -6, 4, 10]
        # both bounds non-periodic
        @test Subzero.in_bounds.(x, y, open_bound, open_bound) .==
            [false, true, false, false, false, true, true, false, false, false]
        # y bounds periodic only
        @test Subzero.in_bounds.(x, y, open_bound, periodic_bound) .==
            [false, true, true, true, true, true, true, true, false, false]
        # x bounds periodic only
        @test Subzero.in_bounds.(x, y, open_bound, periodic_bound) .==
            [true, true, true, false, false, true, true, true, false, false]
        #all bounds periodic 
        @test Subzero.in_bounds.(x, y, periodic_bound, periodic_bound) .==
            [true, true, true, true, true, true, true, true, true, true]

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
        cell_poly = LG.Polygon(cell)
        @test LG.area(cell_poly)::Float64 == 8
        @test LG.GeoInterface.coordinates(cell_poly) == 
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

        # Test aggregate_grid_stress!
        ocean = Subzero.Ocean(grid, 0, 0, 0)
        floe1 = Subzero.Floe([[[1.0, 2], [1, 6], [3, 6], [3, 2], [1, 2]]], 0.5, 0.0)
        # floe is only in cell (7,4) so others will not contribute due to lack of area
        Subzero.aggregate_grid_force!(
            hcat([6, 7, 6, 7, 6, 7], [4, 4, 3, 3, 4, 4]),
            ones(6),
            2ones(6),
            floe1,
            ocean,
            grid,
            open_bound,
            open_bound,
            Threads.SpinLock()
        )
        @test ocean.fx[4, 6] == ocean.fx[3, 6] == ocean.fx[3, 7] == 0
        @test ocean.fx[4,7] == 8
        @test ocean.fy[4, 6] == ocean.fy[3, 6] == ocean.fy[3, 7] == 0
        @test ocean.fy[4, 7] == 16
        @test ocean.si_area[4, 6] == ocean.si_area[3, 6] == ocean.si_area[3, 7] == 0
        @test ocean.si_area[4, 7] == 8 

        floe2 = Subzero.Floe([[[2.0, -4], [2, 0], [6, 0], [6, -4], [2, -4]]], 0.5, 0.0)
        Subzero.aggregate_grid_force!(
            hcat([7, 7, 8, 8, 9, 9], [2, 3, 3, 3, 2, 2]),
            ones(6),
            2ones(6),
            floe2,
            ocean,
            grid,
            periodic_bound,
            periodic_bound,
            Threads.SpinLock()
        )
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
        Subzero.aggregate_grid_force!(
            hcat([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6]),
            ones(6),
            2ones(6),
            floe3,
            ocean,
            grid,
            open_bound,
            open_bound,
            Threads.SpinLock()
        )
        @test ocean.fx[5, 10] == 4
        @test ocean.fx[5, 11] == 2
        @test sum(ocean.fx) == 6

        ocean = Subzero.Ocean(grid, 0, 0, 0)
        Subzero.aggregate_grid_force!(
            hcat([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6]),
            ones(6),
            2ones(6),
            floe3,
            ocean,
            grid,
            periodic_bound,
            open_bound,
            Threads.SpinLock()
        )

        @test ocean.fx[1, 10]  == 8
        @test ocean.fx[2, 11]  == 2
        @test ocean.fx[2, 10] == ocean.fx[1, 11] == 4
        @test sum(ocean.fx) == 18

        ocean = Subzero.Ocean(grid, 0, 0, 0)
        Subzero.aggregate_grid_force!(
            hcat([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6]),
            ones(6),
            2ones(6),
            floe3,
            ocean,
            grid,
            open_bound,
            periodic_bound,
            Threads.SpinLock()
        )
        @test ocean.fx[5, 10] == ocean.fx[5, 1] == 4
        @test ocean.fx[5, 2] == 2
        @test sum(ocean.fx) == 10

        ocean = Subzero.Ocean(grid, 0, 0, 0)
        Subzero.aggregate_grid_force!(
            hcat([10, 10, 11, 11, 12, 12], [5, 6, 5, 6, 5, 6]),
            ones(6),
            2ones(6),
            floe3,
            ocean,
            grid,
            periodic_bound,
            periodic_bound,
            Threads.SpinLock()
        )
        @test ocean.fx[1, 10] == ocean.fx[1, 1] == 8
        @test ocean.fx[2, 10] == ocean.fx[2, 1] == ocean.fx[1, 2] == 4
        @test ocean.fx[2, 2] == 2
        @test sum(ocean.fx) == 30
    end
    @testset "OA Forcings" begin
        # set up model and floe
        grid = Subzero.RegRectilinearGrid(
            Float64,
            (-1e5, 1e5),
            (-1e5, 1e5),
            1e4,
            1e4,
        )
        zonal_ocean = Subzero.Ocean(grid, 1.0, 0.0, 0.0)
        zero_atmos = Subzero.Atmos(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), fill(-20.0, grid.dims .+ 1))
        domain = Subzero.Domain(Subzero.CollisionBoundary(grid, Subzero.North()), Subzero.CollisionBoundary(grid, Subzero.South()),
                        Subzero.CollisionBoundary(grid, Subzero.East()), Subzero.CollisionBoundary(grid, Subzero.West()))
        floe = Subzero.Floe([[[-1.75e4, 5e4], [-1.75e4, 7e4], [-1.25e4, 7e4], 
        [-1.25e4, 5e4], [-1.75e4, 5e4]]], 0.25, 0.0)
        area = floe.area
        # Standard Monte Carlo points for below floe - used to compare with MATLAB
        jldopen("inputs/test_mc_points.jld2", "r") do f
            floe.mc_x = f["X"]
            floe.mc_y = f["Y"]
        end
        
        modulus = 1.5e3*(sqrt(area) + sqrt(area))
        consts = Constants(E = modulus)
        spinlock = Threads.SpinLock()

        # stationary floe, uniform zonal ocean flow
        model1 = Model(
            grid,
            zonal_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model1,
            consts,
            CouplingSettings(Δd = 2),
            spinlock,
        )
        @test isapprox(model1.floes[1].fxOA/area, 2.9760, atol = 1e-3)
        @test isapprox(model1.floes[1].fyOA/area, 0.8296, atol = 1e-3)
        @test isapprox(model1.floes[1].trqOA/area, -523.9212, atol = 1e-3)

        # stationary floe, uniform meridional ocean flow
        meridional_ocean = Subzero.Ocean(grid, 0.0, 1.0, 0.0)
        model2 = Subzero.Model(
            grid,
            meridional_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model2,
            consts,
            CouplingSettings(Δd = 2),
            spinlock,
        )
        @test isapprox(model2.floes[1].fxOA/area, -0.8296, atol = 1e-3)
        @test isapprox(model2.floes[1].fyOA/area, 2.9760, atol = 1e-3)
        @test isapprox(model2.floes[1].trqOA/area, 239.3141, atol = 1e-3)

        # moving floe, uniform 0 ocean flow
        zero_ocean = Subzero.Ocean(grid, 0.0, 0.0, 0.0)
        floe3 = deepcopy(floe)
        floe3.u = 0.25
        floe3.v = 0.1
        model3 = Subzero.Model(
            grid,
            zero_ocean,
            zero_atmos,
            domain,
            StructArray([floe3]),
        )
        Subzero.timestep_coupling!(
            model3,
            consts,
            CouplingSettings(Δd = 2),
            spinlock,
        )
        @test isapprox(model3.floes[1].fxOA/area, -0.1756, atol = 1e-3)
        @test isapprox(model3.floes[1].fyOA/area, -0.1419, atol = 1e-3)
        @test isapprox(model3.floes[1].trqOA/area, 29.0465, atol = 1e-3)

        # rotating floe, uniform 0 ocean flow
        floe4 = deepcopy(floe)
        floe4.ξ = 0.05
        model4 = Subzero.Model(
            grid,
            zero_ocean,
            zero_atmos,
            domain,
            StructArray([floe4]),
        )
        Subzero.timestep_coupling!(
            model4,
            consts,
            CouplingSettings(Δd = 2),
            spinlock,
        )
        @test isapprox(model4.floes[1].fxOA/area, 1.91887860e4, atol = 1e-3)
        @test isapprox(model4.floes[1].fyOA/area, 5.9577026e3, atol = 1e-3)
        @test isapprox(model4.floes[1].trqOA/area, -1.9773119198332e9, atol = 1e-3)
        
        # stationary floe, diagonal atmos flow
        diagonal_atmos = Subzero.Atmos(grid, -1, -0.5, 0.0)
        model5 = Subzero.Model(
            grid,
            zero_ocean,
            diagonal_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model5,
            consts,
            CouplingSettings(Δd = 2),
            spinlock,
        )
        @test isapprox(model5.floes[1].fxOA/area, -0.0013, atol = 1e-3)
        @test isapprox(model5.floes[1].fyOA/area, -6.7082e-4, atol = 1e-3)
        @test isapprox(model5.floes[1].trqOA/area, 0.2276, atol = 1e-3)

        # non-uniform ocean flow, zero atmos, stationary floe
        xgrid, ygrid = Subzero.grids_from_lines(grid.xg, grid.yg)
        psi_ocn = 0.5e4*(sin.(4*(π/4e5).*xgrid) .* sin.(4*(π/4e5).*ygrid))
        non_unif_uocn = zeros(size(xgrid))
        non_unif_uocn[2:end, :] = -1e-4*(psi_ocn[2:end, :] .- psi_ocn[1:end-1, :])
        non_unif_vocn = zeros(size(ygrid))
        non_unif_vocn[:, 2:end] =  1e-4*(psi_ocn[:, 2:end] .- psi_ocn[:, 1:end-1])
        non_unif_ocean = Subzero.Ocean(non_unif_uocn, non_unif_vocn, zeros(size(xgrid)))
        model6 = Subzero.Model(
            grid,
            non_unif_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model6,
            consts,
            CouplingSettings(),
            spinlock,
        )
        @test isapprox(model6.floes[1].fxOA/area, -0.0182, atol = 1e-3)
        @test isapprox(model6.floes[1].fyOA/area, 0.0392, atol = 1e-3)
        @test isapprox(model6.floes[1].trqOA/area, 23.6399, atol = 1e-3)

        # non-uniform atmos flow, zero ocean, stationary floe
        non_unif_atmos = Subzero.Atmos(non_unif_uocn, non_unif_vocn, zeros(size(xgrid)))
        model7 = Subzero.Model(
            grid,
            zero_ocean,
            non_unif_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model7,
            consts,
            CouplingSettings(),
            spinlock,
        )
        @test isapprox(model7.floes[1].fxOA/area, -1.5378e-6, atol = 1e-8)
        @test isapprox(model7.floes[1].fyOA/area, 1.61516e-5, atol = 1e-7)
        @test isapprox(model7.floes[1].trqOA/area, 7.528529e-4, atol = 1e-6)

        # moving floe, non-uniform ocean, non-uniform atmos
        floe8 = deepcopy(floe)
        floe8.u = 0.5
        floe8.v = -0.5
        model8 = Subzero.Model(
            grid,
            non_unif_ocean,
            non_unif_atmos,
            domain,
            StructArray([floe8]),
        )
        Subzero.timestep_coupling!(
            model8,
            consts,
            CouplingSettings(),
            spinlock,
        )
        @test isapprox(model8.floes[1].fxOA/area, -1.6300, atol = 1e-3)
        @test isapprox(model8.floes[1].fyOA/area, 1.1240, atol = 1e-3)
        @test isapprox(model8.floes[1].trqOA/area, 523.2361, atol = 1e-3)
    end
end