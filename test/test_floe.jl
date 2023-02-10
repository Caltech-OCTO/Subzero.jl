@testset "Floe" begin
    # test generating monte carlo points
    file = jldopen("inputs/floe_shapes.jld2", "r")
    floe1_coords = file["floe_vertices"][1]
    poly1 = LibGEOS.Polygon(Subzero.valid_polyvec!(floe1_coords))
    centroid1 = LibGEOS.GeoInterface.coordinates(LibGEOS.centroid(poly1))
    xo, yo = Subzero.seperate_xy(Subzero.translate(floe1_coords, -centroid1))
    rmax = sqrt(maximum([sum(xo[i]^2 + yo[i]^2) for i in eachindex(xo)]))
    area = LibGEOS.area(poly1)
    mc_x, mc_y, alive = Subzero.generate_mc_points(1000, xo, yo, rmax, area, true, Xoshiro(1))
    @test length(mc_x) == length(mc_y) && length(mc_x) > 0
    in_on = inpoly2(hcat(mc_x, mc_y), hcat(xo, yo))
    mc_in = in_on[:, 1] .|  in_on[:, 2]
    @test all(mc_in)
    @test abs(sum(mc_in)/1000 * 4 * rmax^2 - area)/area < 0.1
    @test alive
    # Test that random number generator is working
    mc_x2, mc_y2, alive2 = Subzero.generate_mc_points(1000, xo, yo, rmax, area, true, Xoshiro(1))
    @test all(mc_x .== mc_x2)
    @test all(mc_y .== mc_y2)
    @test alive == alive2

    # Test InteractionFields enum
    interactions = range(1, 7)'
    @test interactions[floeidx] == 1
    @test interactions[xpoint] == 4
    @test interactions[overlap] == 7

    # Test with coords inputs
    floe_from_coords = Floe(floe1_coords, 0.5, 0.01, u = 0.2, rng = Xoshiro(1))
    @test typeof(floe_from_coords) <: Floe
    @test floe_from_coords.u == 0.2
    @test 0.49 <= floe_from_coords.height <= 0.51
    @test floe_from_coords.centroid == centroid1
    @test floe_from_coords.area == area
    @test floe_from_coords.alive
    
    # Test with polygon input
    floe_from_poly = Floe(poly1, 0.5, 0.01, v = -0.2, rng = Xoshiro(1))
    @test typeof(floe_from_poly) <: Floe
    @test floe_from_poly.u == 0.0
    @test floe_from_poly.v == -0.2
    @test 0.49 <= floe_from_poly.height <= 0.51
    @test floe_from_poly.centroid == centroid1
    @test floe_from_poly.area == area
    @test floe_from_poly.alive

    # Test poly_to_floes

    # Test initialize_floe_field from coord list

    # Test initialize_floe_field with voronoi
end