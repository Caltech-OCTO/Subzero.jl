@testset "Floe Operations" begin
    exterior = [[0.0, 1.0],  [0.0, 0.0],  [1.0, 0.0],  [1.0, 1.0], [0.0, 1.0]]
    hole = [[0.2, 0.7], [0.2, 0.2], [0.7, 0.2], [0.7, 0.7], [0.2, 0.7]]
    poly_nohole = LibGEOS.Polygon([exterior])
    poly_hole = LibGEOS.Polygon([exterior, hole])
    @test !Subzero.hashole(poly_nohole)
    @test Subzero.hashole(poly_hole)
end