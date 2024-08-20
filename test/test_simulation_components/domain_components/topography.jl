@testset "Topography" begin
    coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]]
    poly = Subzero.make_polygon(coords)
    # Polygon Constructor
    topo1 = Subzero.TopographyElement(; poly)
    # @test topo1.coords == coords
    @test topo1.centroid == [0.5, 0.5]
    @test topo1.rmax == sqrt(0.5)
    topo32 = Subzero.TopographyElement(Float32; poly)
    @test typeof(topo32) == Subzero.TopographyElement{Float32}
    # @test typeof(topo32.coords) == Subzero.PolyVec{Float32}
    # Coords Constructor
    topo2 = Subzero.TopographyElement(Float64; coords)
    # @test topo2.coords == coords
    @test topo2.centroid == [0.5, 0.5]
    @test topo2.rmax == sqrt(0.5)
    # Basic constructor
    topo3 = TopographyElement{Float64}(
        Subzero.make_polygon(coords),
        [0.5, 0.5],
        sqrt(0.5),
    )
    # @test topo3.coords == coords
    # check when radius is less than  or equal to 0
    @test_throws ArgumentError TopographyElement{Float64}(
        Subzero.make_polygon(coords),
        [0.5, 0.5],
        -sqrt(0.5),
    )

    # Create field of topography
    coords_w_hole = [
        [[0.5, 10.0], [0.5, 0.0], [10.0, 0.0], [10.0, 10.0], [0.5, 10.0]],
        [[2.0, 8.0], [2.0, 4.0], [8.0, 4.0], [8.0, 8.0], [2.0, 8.0]]
        ]
    topo_field_64 = initialize_topography_field(Float64; coords = [coords, coords_w_hole])
    @test length(topo_field_64) == 2
    @test typeof(topo_field_64) <: StructArray{TopographyElement{Float64}}
    # @test !Subzero.hashole(topo_field_64.coords[2])

    topo_field_32 = initialize_topography_field(Float32; coords = [coords, coords_w_hole])
    @test typeof(topo_field_32) <: StructArray{TopographyElement{Float32}}
end