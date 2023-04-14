@testset "Simplification" begin
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
            CouplingSettings(),
            2,
            Xoshiro(1),
        )
        @test max_id == 2
        @test f1.coords == coords1
        @test f2.coords == coords2

        # Test two floes intersecting -> will fuse into floe1 since same size
        Subzero.translate!(coords2, -13.0, 0.0)
        f2 = Floe(coords2, 0.5, 0.0)
        f1.id = 1
        f2.id = 2
        mass_tot = f1.mass + f2.mass
        f1.u = 0.1
        f1.v = 0.1
        f1.ξ = 0.1
        f1.p_dudt = 0.01
        f2.u = 0.1
        f2.v = 0.2
        f2.ξ = 0
        f2.p_dudt = 0.02
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

        max_id = Subzero.fuse_floes!(
            f1,
            f2,
            Constants(),
            CouplingSettings(),
            2,
            Xoshiro(1),
        )
        @test max_id == 3
        @test f1.status.tag == Subzero.active
        @test f2.status.tag == Subzero.remove
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
        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-6)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-6)
        @test isapprox(
            spin_momentum_init + angular_momentum_init,
            spin_momentum_after + angular_momentum_after,
            atol = 1e-6,
        )
        # Test two floes intersecting -> will fuse into floe2 since bigger
    end
end