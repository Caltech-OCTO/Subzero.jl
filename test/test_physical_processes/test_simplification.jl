@testset "Simplification" begin
    @testset "Dissolve Floes" begin
        grid = RegRectilinearGrid(
            Float64,
            (-1e5, 1e5),
            (0.0, 1e5),
            1e4,
            1e4,
        )
        domain = Subzero.Domain(
            CollisionBoundary(grid, North()),
            CollisionBoundary(grid, South()),
            PeriodicBoundary(grid, East()),
            PeriodicBoundary(grid, West()),
        )
        height = 0.25
        ρi = 920.0
        coords = [[
            [0.0, 5e4],
            [0.0, 8e4],
            [3e4, 8e4],
            [3e4, 5e4],
            [0.0, 5e4],
        ]]
        mass = 9e8 * height * ρi
        dissolve = zeros(Float64, 10, 20)  # 20x20 ocean grid
        # Add 2 floes in the middle of the grid -> masses added to cell
        floe = Floe(coords, height, 0.0,)
        Subzero.dissolve_floe!(floe, grid, domain, dissolve)
        @test dissolve[7, 12] == mass
        @test floe.status.tag == Subzero.remove
        floe = Floe(Subzero.translate(coords, 2.5e3, 2.5e3), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolve)
        @test dissolve[7, 12] == 2mass
        @test floe.status.tag == Subzero.remove
        # Add floe over periodic bound -> mass added to cell wrapped around grid
        floe = Floe(Subzero.translate(coords, 9e4, 0.0), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolve)
        @test dissolve[7, 1] == mass
        @test floe.status.tag == Subzero.remove
        floe = Floe(Subzero.translate(coords, -1.2e5, 0.0), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolve)
        @test dissolve[7, 20] == mass
        @test floe.status.tag == Subzero.remove
        total_mass = sum(dissolve)
        # Add floe over non-periodic bound -> mass not added since out of bounds
        floe = Floe(Subzero.translate(coords, 0.0, 6e4), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolve)
        @test total_mass == sum(dissolve)  # nothing was added
        @test floe.status.tag == Subzero.remove
        floe = Floe(Subzero.translate(coords, 0.0, -7e4), height, 0.0)
        Subzero.dissolve_floe!(floe, grid, domain, dissolve)
        @test total_mass == sum(dissolve)  # nothing was added
        @test floe.status.tag == Subzero.remove
    end

    @testset "Fuse Floes" begin
        # Test two floes not intersecting -> will not fuse
        coords1 = [[
            [0.0, 0.0],
            [0.0, 10.0],
            [10.0, 10.0],
            [10.0, 0.0],
            [0.0, 0.0],
        ]]
        coords2 = deepcopy(coords1)
        Subzero.translate!(coords2, 20.0, 0.0)
        f1 = Floe(coords1, 0.5, 0.0)
        f2 = Floe(coords2, 0.5, 0.0)
        max_id = Subzero.fuse_floes!(
            f1,
            f2,
            Subzero.Constants(),
            10,
            CouplingSettings(),
            2,
            Xoshiro(1),
        )
        @test max_id == 2
        @test f1.coords == coords1
        @test f2.coords == coords2

        # Test two floes intersecting -> will fuse into floe1 since same size
        Subzero.Subzero.translate!(coords2, -13.0, 0.0)
        f2 = Floe(coords2, 0.75, 0.0)
        f1.id = 1
        f2.id = 2
        mass_tot = f1.mass + f2.mass
        f1.u = 0.1
        f1.v = 0.1
        f1.ξ = 0.1
        f1.p_dxdt = 0.002
        f1.p_dydt = 0.08
        f1.p_dαdt = 0.04
        f1.p_dudt = 0.01
        f2.p_dvdt = -0.007
        f2.p_dξdt = 0.01
        f2.u = 0.1
        f2.v = 0.2
        f2.ξ = 0
        f2.p_dudt = 0.02
        f2.p_dvdt = -0.005
        f2.p_dξdt = 0.05
        stress1_init = f1.stress
        x_momentum_init, y_momentum_init = Subzero.calc_linear_momentum(
            [f1.u, f2.u],
            [f1.v, f2.v],
            [f1.mass, f2.mass],
        )
        spin_momentum_init, angular_momentum_init = Subzero.calc_angular_momentum(
            [f1.u, f2.u],
            [f1.v, f2.v],
            [f1.mass, f2.mass],
            [f1.ξ, f2.ξ],
            [f1.moment, f2.moment],
            [f1.centroid[1], f2.centroid[1]],
            [f1.centroid[2], f2.centroid[2]],
        )
        p_x_momentum_init, p_y_momentum_init = Subzero.calc_linear_momentum(
            [f1.p_dxdt, f2.p_dxdt],
            [f1.p_dydt, f2.p_dydt],
            [f1.mass, f2.mass],
        )
        p_spin_momentum_init, p_angular_momentum_init = Subzero.calc_angular_momentum(
            [f1.p_dxdt, f2.p_dxdt],
            [f1.p_dydt, f2.p_dydt],
            [f1.mass, f2.mass],
            [f1.p_dαdt, f2.p_dαdt],
            [f1.moment, f2.moment],
            [f1.centroid[1] - 10 * f1.p_dxdt, f2.centroid[1] - 10 * f2.p_dxdt],
            [f1.centroid[2] - 10 * f1.p_dydt, f2.centroid[2] - 10 * f2.p_dydt],
        )

        max_id = Subzero.fuse_floes!(
            f1,
            f2,
            Constants(),
            10,
            CouplingSettings(),
            2,
            Xoshiro(1),
        )
        @test max_id == 3
        @test f1.status.tag == Subzero.active
        @test f2.status.tag == Subzero.Subzero.remove
        @test f1.area == 170
        @test f1.mass == mass_tot  # conservation of mass
        @test f1.parent_ids == [1, 2]
        @test f1.id == 3
        # conservation of momentum
        x_momentum_after, y_momentum_after = Subzero.calc_linear_momentum(
            [f1.u],
            [f1.v],
            [f1.mass],
        )
        spin_momentum_after, angular_momentum_after = Subzero.calc_angular_momentum(
            [f1.u],
            [f1.v],
            [f1.mass],
            [f1.ξ],
            [f1.moment],
            [f1.centroid[1]],
            [f1.centroid[2]],
        )
        p_x_momentum_after, p_y_momentum_after = Subzero.calc_linear_momentum(
            [f1.p_dxdt],
            [f1.p_dydt],
            [f1.mass],
        )
        p_spin_momentum_after, p_angular_momentum_after = Subzero.calc_angular_momentum(
            [f1.p_dxdt],
            [f1.p_dydt],
            [f1.mass],
            [f1.p_dαdt],
            [f1.moment],
            [f1.centroid[1] - 10 * f1.p_dxdt],
            [f1.centroid[2] - 10 * f1.p_dydt],
        )
        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-12)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-12)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-12)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-12)
        @test isapprox(
            spin_momentum_init + angular_momentum_init,
            spin_momentum_after + angular_momentum_after,
            atol = 1e-6,
        )
        # this is much larger due to all of the uncertianty with centroids
        @test isapprox(
            p_spin_momentum_init + p_angular_momentum_init,
            p_spin_momentum_after + p_angular_momentum_after,
            atol = 1e3,
        )
        @test mean(f1.stress_history.cb) == f1.stress_history.total/1000 == f1.stress
        @test f1.stress == (stress1_init * (f2.mass - mass_tot) .+ f2.stress * f2.mass) / mass_tot

        # Test two floes intersecting -> will fuse into floe2 since bigger
        f1 = Floe(coords1, 0.5, 0.0)
        f3 = Floe(
            [[
                [0.0, 0.0],
                [0.0, 20.0],
                [20.0, 20.0],
                [20.0, 0.0],
                [0.0, 0.0],
            ]],
            0.55,
            0.0,
        )
        max_id = Subzero.fuse_floes!(
            f1,
            f3,
            Constants(),
            10,
            CouplingSettings(),
            3,
            Xoshiro(1),
        )
        f1.id = 1
        f2.id = 2
        @test f1.status.tag == Subzero.Subzero.remove
        @test f3.status.tag == Subzero.active
    end
end