@testset "Floe" begin
    # test generating monte carlo points
    file = jldopen("Subzero/examples/floe_shapes.jld2", "r")
    floe1_coords = file["floe_vertices"][1]
    poly1 = LibGEOS.Polygon(floe1_coords)
    centroid1 = LibGEOS.GeoInterface.coordinates(LibGEOS.centroid(poly1))
    xo, yo = Sybzero.seperate_xy(Subzero.translate(floe1_coords, -centroid))
    rmax = sqrt(maximum([sum(xo[i]^2 + yo[i]^2) for i in eachindex(xo)]))
    area = LibGEOS.area(poly1)
    mc_x, mc_y, alive = generate_mc_points(1000, xo, yo, rmax, area, true, Xoshiro(1))
    @test length(mc_x) == length(mc_y) && length(mc_x) > 0
    in_on = inpoly2(hcat(mc_x, mc_y), hcat(xo, yo))
    mc_in = in_on[:, 1] .|  in_on[:, 2]
    @test all(mc_in)
    @test abs(sum(mc_in)/1000 * 4 * rmax^2 - area)/area


    # test with polygon inputs
    # test with coords input
end