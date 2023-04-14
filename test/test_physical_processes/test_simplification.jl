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
        translate!(coords2, 10.0, 0.0)
        
    end
end