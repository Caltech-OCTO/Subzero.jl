@testset "Simulation Creation" begin
    @testset "Floe Area Ratio" begin
        @test Subzero.floe_grid_bounds(collect(0:2:10), 5, 1.25) == (2, 5)
        @test Subzero.floe_grid_bounds(collect(0:2:10), 5, 2.25) == (2, 5)
        @test Subzero.floe_grid_bounds(collect(0:2:10), 0.75, 1.25) == (1, 2)
        @test Subzero.floe_grid_bounds(collect(0:2:10), 9.25, 1.25) == (5, 6)

        cell_poly = LibGEOS.Polygon([[[0., 10.], [0., 0.], [10., 0.],
                                      [10., 10.], [0., 10.]]])
        floe_poly1 = LibGEOS.Polygon([[[2., 8.], [2., 2.], [8., 2.],
                                         [8., 8.], [2., 8.]]])
        floe_poly2 = LibGEOS.Polygon([[[9., 5.], [9., 2.], [12., 2.],
                                         [12., 5.], [9., 5.]]])
        @test Subzero.cell_area_ratio(cell_poly, floe_poly1) == 0.36
        @test Subzero.cell_area_ratio(cell_poly, floe_poly2) == 0.03
        @test Subzero.cell_area_ratio(cell_poly,
                Subzero.translate(floe_poly2, [1.0, 0.0])) == 0.0
        floe_poly3 = LibGEOS.Polygon([[[3., 8.], [3., 3.], [8., 3.],
                                       [8., 8.],[3., 8.]]])
        floe3 = Subzero.Floe(floe_poly3, 0.25, 0.0)
        area_ratio3, _, _, idx1 = Subzero.floe_area_ratio(floe3, collect(-5.:5.:10.), collect(0.:5.:10.))
        @test area_ratio3[1] == 4/25
        @test area_ratio3[2] == 6/25
        @test area_ratio3[3] == 6/25
        @test area_ratio3[4] == 9/25

        floe_poly4 = LibGEOS.Polygon([[[2., 2.], [8., 2.], [8., 8.], [2., 2.]]])
        floe4 = Subzero.Floe(floe_poly4, 0.25, 0.0)
        area_ratio4, _, _, idx1 = Subzero.floe_area_ratio(floe4, collect(-5.:5.:10.), collect(0.:5.:10.))
        @test area_ratio4[1] == 0.18
        @test area_ratio4[2] == 0.36
        @test area_ratio4[3] == 0.18
    end

    @testset "Collisions" begin
        coords1 = [[[0.0, 0.0], [0.0, 10000.0], [5000, 10000], [5000, 0], [0, 0]]]
        coords2 = [[[4998.0, 4000.0], [4998.0, 8000.0], [10000, 8000], [10000, 4000], [4998.0, 4000]]]
        floe1 = Subzero.Floe(LibGEOS.Polygon(coords1), 0.25, 0.0)
        floe2 = Subzero.Floe(LibGEOS.Polygon(coords2), 0.25, 0.0)

        coords1 = [[[0., -1000], [0, 0], [1000, 0], [0, -1000]]]
        coords2 = [[[925., -700], [900, 10], [950, 20], [975, -600], [925, -700]]]
    end

end