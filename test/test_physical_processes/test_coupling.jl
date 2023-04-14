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
        open_open_answers = [false, true, true, false, false, true, true, true, false, false]
        periodic_open_answers = [false, true, true, true, true, true, true, true, false, false]
        open_periodic_answers = [true, true, true, false, false, true, true, true, true, false]
        periodic_periodic_answers = [true, true, true, true, true, true, true, true, true, true]
        # both bounds non-periodic
        for i in eachindex(x)
            @test Subzero.in_bounds(x[i], y[i], grid, open_bound, open_bound) == open_open_answers[i]
            @test Subzero.in_bounds(x[i], y[i], grid, open_bound, periodic_bound) == open_periodic_answers[i]
            @test Subzero.in_bounds(x[i], y[i], grid, periodic_bound, open_bound) == periodic_open_answers[i]
            @test Subzero.in_bounds(x[i], y[i], grid, periodic_bound, periodic_bound) == periodic_periodic_answers[i]
        end

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
        function test_floe_to_grid(
            floeidx,
            rows,
            cols,
            τxs,
            τys,
            grid,
            bound1,
            bound2,
            scells,
            settings,
            occupied_cells,
            Δxs,
            Δys,
            sum_τx,
            sum_τy,
            npoints,
        )
            for i in eachindex(rows)
                Subzero.floe_to_grid_info!(
                    floeidx,
                    rows[i],
                    cols[i],
                    τxs[i],
                    τys[i],
                    grid,
                    bound1,
                    bound2,
                    scells,
                    settings,
                )
            end

            for i in eachindex(occupied_cells)
                @test grid.floe_locations[occupied_cells[i]].floeidx[end] == floeidx
                @test grid.floe_locations[occupied_cells[i]].Δx[end] == Δxs[i]
                @test grid.floe_locations[occupied_cells[i]].Δy[end] == Δys[i]
                @test scells[occupied_cells[i]].τx[end] == sum_τx[i]
                @test scells[occupied_cells[i]].τy[end] == sum_τy[i]
                @test scells[occupied_cells[i]].npoints[end] == npoints[i]
            end

            empty_check = true
            for cidx in CartesianIndices(scells)
                if !(cidx in occupied_cells)
                    empty_check = empty_check &&
                        isempty(grid.floe_locations[cidx].floeidx) &&
                        isempty(grid.floe_locations[cidx].Δx) &&
                        isempty(grid.floe_locations[cidx].Δy)
                    empty_check = empty_check &&
                        isempty(scells[cidx].τx) &&
                        isempty(scells[cidx].τy) &&
                        isempty(scells[cidx].npoints)
                end
            end
            @test empty_check
        end
        # Open bounds, everything within bounds
        test_floe_to_grid(
            1,
            [4, 4, 3, 3, 4, 4],
            [7, 7, 6, 6, 7, 7],
            ones(Float64, 6),
            2ones(Float64, 6),
            deepcopy(grid),
            open_bound,
            open_bound,
            Subzero.Ocean(Float64, grid, 0, 0, 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [CartesianIndex(4, 7), CartesianIndex(3, 6)],
            [0.0, 0.0],
            [0.0, 0.0],
            [-4, -2],
            [-8, -4],
            [4, 2],
        )
        # Periodic bounds, everything within bounds
        test_floe_to_grid(
            2,
            [2, 3, 3, 3, 2, 2],
            [7, 7, 8, 8, 9, 9],
            ones(Float64, 6),
            2ones(Float64, 6),
            deepcopy(grid),
            periodic_bound,
            periodic_bound,
            Subzero.Ocean(Float64, grid, 0, 0, 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(2, 7),
                CartesianIndex(3, 7),
                CartesianIndex(2, 9),
                CartesianIndex(3, 8),
            ],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [-1, -1, -2, -2],
            [-2, -2, -4, -4],
            [1, 1, 2, 2],
        )
        # Periodic bounds in NS, floe out of bounds - moved to other side of grid
        test_floe_to_grid(
            3,
            [4, 5, 6, 5, 6, 5, 6],
            [10, 10, 10, 11, 11, 11, 11],
            ones(Float64, 7),
            2ones(Float64, 7),
            deepcopy(grid),
            periodic_bound,
            open_bound,
            Subzero.Ocean(Float64, grid, 0, 0, 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(1, 10),
                CartesianIndex(1, 11),
                CartesianIndex(2, 10),
                CartesianIndex(2, 11),
                CartesianIndex(4, 10),
            ],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [-16.0, -16.0, -16.0, -16.0, 0.0],
            [-1, -2, -1, -2, -1],
            [-2, -4, -2, -4, -2],
            [1, 2, 1, 2, 1],
        )
        # Periodic bounds in EW, floe out of bounds - moved to other side of grid
        test_floe_to_grid(
            4,
            [4, 5, 5, 5, 4],
            [11, 11, 12, 12, 11],
            ones(Float64, 5),
            2ones(Float64, 5),
            deepcopy(grid),
            open_bound,
            periodic_bound,
            Subzero.Ocean(Float64, grid, 0, 0, 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(4, 1),
                CartesianIndex(5, 1),
                CartesianIndex(5, 2),
            ],
            [-20.0, -20.0, -20.0],
            [0.0, 0.0, 0.0],
            [-2, -1, -2],
            [-4, -2, -4],
            [2, 1, 2],
        )
        # Periodic bounds in both direction, floe out of bounds
        test_floe_to_grid(
            2,
            [0, -1, -2, 1, -1],
            [0, -1, -1, 1, -1],
            -ones(Float64, 5),
            -2ones(Float64, 5),
            deepcopy(grid),
            periodic_bound,
            periodic_bound,
            Subzero.Ocean(Float64, grid, 0, 0, 0).scells,
            CouplingSettings(two_way_coupling_on = true),
            [
                CartesianIndex(1, 1),
                CartesianIndex(4, 10),
                CartesianIndex(3, 9),
                CartesianIndex(2, 9),
            ],
            [0.0, 20.0, 20.0, 20.0],
            [0.0, 16.0, 16.0, 16.0],
            [1, 1, 2, 1],
            [2, 2, 4, 2],
            [1, 1, 2, 1],
        )
    end
    @testset "OA Forcings" begin
    #     # set up model and floe
        grid = Subzero.RegRectilinearGrid(
            Float64,
            (-1e5, 1e5),
            (-1e5, 1e5),
            1e4,
            1e4,
        )
        zonal_ocean = Subzero.Ocean(Float64, grid, 1.0, 0.0, 0.0)
        zero_atmos = Subzero.Atmos(grid, 0.0, 0.0, -20.0)
        domain = Subzero.Domain(
            Subzero.CollisionBoundary(grid, Subzero.North()),
            Subzero.CollisionBoundary(grid, Subzero.South()),
            Subzero.CollisionBoundary(grid, Subzero.East()),
            Subzero.CollisionBoundary(grid, Subzero.West()),
        )
        floe = Subzero.Floe(
            [[
                [-1.75e4, 5e4],
                [-1.75e4, 7e4],
                [-1.25e4, 7e4], 
                [-1.25e4, 5e4],
                [-1.75e4, 5e4],
            ]],
            0.25,
            0.0,
        )
        area = floe.area
    #     # Standard Monte Carlo points for below floe - used to compare with MATLAB
        jldopen("inputs/test_mc_points.jld2", "r") do f
            floe.mc_x = f["X"]
            floe.mc_y = f["Y"]
        end
        
        modulus = 1.5e3*(sqrt(area) + sqrt(area))
        consts = Constants(E = modulus)

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
            10,
            consts,
            CouplingSettings(Δd = 2),
        )
        @test isapprox(model1.floes[1].fxOA/area, 2.9760, atol = 1e-3)
        @test isapprox(model1.floes[1].fyOA/area, 0.8296, atol = 1e-3)
        @test isapprox(model1.floes[1].trqOA/area, -523.9212, atol = 1e-3)

    # stationary floe, uniform meridional ocean flow
        meridional_ocean = Subzero.Ocean(Float64, grid, 0.0, 1.0, 0.0)
        model2 = Subzero.Model(
            grid,
            meridional_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model2,
            10,
            consts,
            CouplingSettings(Δd = 2),
        )
        @test isapprox(model2.floes[1].fxOA/area, -0.8296, atol = 1e-3)
        @test isapprox(model2.floes[1].fyOA/area, 2.9760, atol = 1e-3)
        @test isapprox(model2.floes[1].trqOA/area, 239.3141, atol = 1e-3)

    # moving floe, uniform 0 ocean flow
        zero_ocean = Subzero.Ocean(Float64, grid, 0.0, 0.0, 0.0)
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
            10,
            consts,
            CouplingSettings(Δd = 2),
        )
        @test isapprox(model3.floes[1].fxOA/area, -0.1756, atol = 1e-3)
        @test isapprox(model3.floes[1].fyOA/area, -0.1419, atol = 1e-3)
        @test isapprox(model3.floes[1].trqOA/area, 29.0465, atol = 1e-1)
        
        # stationary floe, diagonal atmos flow
        diagonal_atmos = Subzero.Atmos(grid, -1, -0.5, 0.0)
        model4 = Subzero.Model(
            grid,
            zero_ocean,
            diagonal_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model4,
            10,
            consts,
            CouplingSettings(Δd = 2),
        )
        @test isapprox(model4.floes[1].fxOA/area, -0.0013, atol = 1e-3)
        @test isapprox(model4.floes[1].fyOA/area, -6.7082e-4, atol = 1e-3)
        @test isapprox(model4.floes[1].trqOA/area, 0.2276, atol = 1e-3)

        # non-uniform ocean flow, zero atmos, stationary floe
        xgrid, ygrid = Subzero.grids_from_lines(grid.xg, grid.yg)
        psi_ocn = 0.5e4*(sin.(4*(π/4e5).*xgrid) .* sin.(4*(π/4e5).*ygrid))
        non_unif_uocn = zeros(size(xgrid))
        non_unif_uocn[2:end, :] = -1e-4*(psi_ocn[2:end, :] .- psi_ocn[1:end-1, :])
        non_unif_vocn = zeros(size(ygrid))
        non_unif_vocn[:, 2:end] = 1e-4*(psi_ocn[:, 2:end] .- psi_ocn[:, 1:end-1])
        non_unif_ocean = Subzero.Ocean(
            Float64,
            non_unif_uocn,
            non_unif_vocn,
            zeros(size(xgrid)),
        )
        model5 = Subzero.Model(
            grid,
            non_unif_ocean,
            zero_atmos,
            domain,
            StructArray([deepcopy(floe)]),
        )
        Subzero.timestep_coupling!(
            model5,
            10,
            consts,
            CouplingSettings(),
        )
        @test isapprox(model5.floes[1].fxOA/area, -0.0182, atol = 1e-3)
        @test isapprox(model5.floes[1].fyOA/area, 0.0392, atol = 1e-3)
        @test isapprox(model5.floes[1].trqOA/area, 23.6399, atol = 1e-3)


        # moving floe, non-uniform ocean, non-uniform atmos
        non_unif_atmos = Subzero.Atmos(
            non_unif_uocn,
            non_unif_vocn,
            zeros(size(xgrid)),
        )
        floe6 = deepcopy(floe)
        floe6.u = 0.5
        floe6.v = -0.5
        model6 = Subzero.Model(
            grid,
            non_unif_ocean,
            non_unif_atmos,
            domain,
            StructArray([floe6]),
        )
        Subzero.timestep_coupling!(
            model6,
            10,
            consts,
            CouplingSettings(),
        )
        @test isapprox(model6.floes[1].fxOA/area, -1.6300, atol = 1e-3)
        @test isapprox(model6.floes[1].fyOA/area, 1.1240, atol = 1e-3)
        @test isapprox(model6.floes[1].trqOA/area, 523.2361, atol = 2e-1)
    end
end