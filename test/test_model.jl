@testset "Model Creation" begin
    # Grid Creation
    @testset "Grid" begin
        # Large grid default constructor
        xg = collect(-4e5:1e4:4e5)
        yg = collect(-3e5:1e4:3e5)
        xc = collect(-3.95e5:1e4:3.95e5)
        yc = collect(-2.95e5:1e4:2.95e5)
        gbig = Subzero.RegRectilinearGrid((length(yc), length(xc)), xg, yg, xc, yc)
        @test gbig.dims == (60, 80)
        @test gbig.xg == xg
        @test gbig.yg == yg
        @test gbig.xc == xc
        @test gbig.yc == yc
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.RegRectilinearGrid((60, 80), xg, yg, Float64[0.0], yc)
        @test_throws ArgumentError Subzero.RegRectilinearGrid((60, 80), xg, yg, xc, Float64[0.0])
        @test_throws ArgumentError Subzero.RegRectilinearGrid((60, 80), Float64[0.0], yg, xc, yc)
        @test_throws ArgumentError Subzero.RegRectilinearGrid((60, 80), xg, Float64[0.0], xc, yc)
        
        # Non-square grid using custom constructor
        g1 = Subzero.RegRectilinearGrid(-10, 10, -8, 8, 2, 4)
        @test g1.dims == (4, 10)
        @test g1.xg == collect(-10.0:2:10.0)
        @test g1.yg == collect(-8:4:8)
        @test g1.xc == collect(-9.0:2:9.0)
        @test g1.yc == collect(-6.0:4:6.0)
        @test typeof(g1) == Subzero.RegRectilinearGrid{Float64}
        # Uneven grid size creation (grid cut short) using custom constructor
        g2 = Subzero.RegRectilinearGrid(0.0, 10.5, 0.0, 8.0, 2.5, 2)
        @test g2.dims == (4, 4)
        @test g2.xg == collect(0.0:2.5:10.5)
        @test g2.yg == collect(0.0:2.0:8.0)
        @test g2.xc == collect(1.25:2.5:8.75)
        @test g2.yc == collect(1.0:2.0:7.0)
        # Custom constructor Float32
        @test typeof(Subzero.RegRectilinearGrid(0, 10, 0, 8, 2, 2, Float32)) ==
              Subzero.RegRectilinearGrid{Float32}
        # Grid constructor with dimensions
        g3 = Subzero.RegRectilinearGrid(-10, 10, -8, 8, (4, 10))
        @test g3.dims == (4,10)
        @test g3.xg == g1.xg
        @test g3.yc == g1.yc
    end

    @testset "Ocean" begin
        g = Subzero.RegRectilinearGrid(0, 4e5, 0, 3e5, 1e4, 1e4)
        # Large ocean default constructor
        uocn = fill(3.0, g.dims .+ 1)
        vocn = fill(4.0, g.dims .+ 1)
        tempocn = fill(-2.0, g.dims .+ 1)
        fx = fill(0.0, g.dims .+ 1)
        fy = fx
        si_area = fx
        hflx = fx
        ocn = Subzero.Ocean(uocn, vocn, tempocn, hflx, fx, fy, si_area)
        @test ocn.u == uocn
        @test ocn.v == vocn
        @test ocn.temp == tempocn
        @test ocn.fx == fx == ocn.fy == ocn.si_area == ocn.hflx
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Ocean(Matrix{Float64}(undef, 0, 0), vocn, tempocn, hflx, fx, fy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, Matrix{Float64}(undef, 0, 0), tempocn, hflx, fx, fy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, Matrix{Float64}(undef, 0, 0), hflx, fx, fy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, Matrix{Float64}(undef, 0, 0), fx, fy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, hflx, Matrix{Float64}(undef, 0, 0), fy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, hflx, fx, fy, Matrix{Float64}(undef, 0, 0))
        # Custom constructor
        ocn2 = Subzero.Ocean(g, 3.0, 4.0, -2.0)
        @test ocn.u == ocn2.u
        @test ocn.v == ocn2.v
        @test ocn.temp == ocn2.temp
        @test ocn.fx == ocn2.fx
        @test ocn.si_area == ocn2.si_area
        @test ocn.hflx == ocn2.hflx
        # Custom constructor Float32
        @test typeof(Subzero.Ocean(g, 3.0, 4.0, -2.0, Float32)) ==
              Subzero.Ocean{Float32}
    end

    @testset "Wind" begin
        g = Subzero.RegRectilinearGrid(0, 4e5, 0, 3e5, 1e4, 1e4)
        # Large wind default constructor
        uwind = fill(3.0, g.dims .+ 1)
        vwind = fill(4.0, g.dims .+ 1)
        tempwind = fill(-2.0, g.dims .+ 1)
        wind = Subzero.Wind(uwind, vwind, tempwind)
        @test wind.u == uwind
        @test wind.v == vwind
        @test wind.temp == tempwind
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Wind(Matrix{Float64}(undef, 0, 0), vwind, tempwind)
        @test_throws ArgumentError Subzero.Wind(uwind, Matrix{Float64}(undef, 0, 0), tempwind)
        @test_throws ArgumentError Subzero.Wind(uwind, vwind, Matrix{Float64}(undef, 0, 0))
        # Custom constructor
        wind2 = Subzero.Wind(g, 3.0, 4.0, -2.0)
        @test wind.u == wind2.u
        @test wind.v == wind2.v
        @test wind.temp == wind2.temp
        # Custom constructor Float32
        @test typeof(Subzero.Wind(g, 3.0, 4.0, -2.0, Float32)) ==
              Subzero.Wind{Float32}
    end

    @testset "Boundaries" begin
        # Boundary
        g = Subzero.RegRectilinearGrid(0, 4e5, 0, 3e5, 1e4, 1e4)
        b1 = Subzero.PeriodicBoundary(g, Subzero.North())
        b2 = Subzero.OpenBoundary(g, Subzero.East())
        b3 = Subzero.CollisionBoundary(g, Subzero.West())
        b4 = Subzero.PeriodicBoundary(g, Subzero.South())
        b5 = Subzero.CompressionBoundary(g, Subzero.South(), 1.0)
        @test b1.val == 3e5
        @test typeof(b1) == Subzero.PeriodicBoundary{Subzero.North, Float64}
        @test b1.coords == [[[-2e5, 3e5], [-2e5, 4.5e5], [6e5, 4.5e5], [6e5, 3e5], [-2e5, 3e5]]]
        @test b3.val == 0.0
        @test typeof(b3) == Subzero.CollisionBoundary{Subzero.West, Float64}
        @test b3.coords == [[[-2e5, -1.5e5], [-2e5, 4.5e5], [0.0, 4.5e5], [0.0, -1.5e5], [-2e5, -1.5e5]]]
        @test b5.velocity == 1.0

        # Test changing CompressionBoundary values

        # Test changing CollisionBoundary values

        # Need to add in ability to run with Float32
        #b32 = Subzero.OpenBoundary(Subzero.OpenBC(), 1, Float32)
        #@test typeof(b32.val) == Float32

        # Periodic Compat
        @test !Subzero.periodic_compat(b1, b2)
        @test !Subzero.periodic_compat(b2, b1)
        @test Subzero.periodic_compat(b1, b4)
        @test Subzero.periodic_compat(b2, b3)

        # Test Boundary Coords

        

        # RectangularDomain default constructor
        rdomain1 = Subzero.Domain(b1, b4, b2, b3)
        @test rdomain1.north == b1

        # Domain fails
        # Directions are not correct
        @test_throws MethodError Subzero.Domain(b4, b2, b2, b3)
        @test_throws ArgumentError Subzero.Domain(b1, Subzero.OpenBoundary(g, Subzero.South()), b2, b3)
        #@test_throws MethodError Subzero.Domain(b4, b1, b3, b32)
    end

    @testset "Topography" begin
        coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]]
        poly = LibGEOS.Polygon(coords)
        # Polygon Constructor
        topo1 = Subzero.TopographyElement(poly)
        @test topo1.coords == coords
        topo32 = Subzero.TopographyElement(poly, Float32)
        @test typeof(topo32) == Subzero.TopographyElement{Float32}
        @test typeof(topo32.coords) == Subzero.PolyVec{Float32}
        # Coords Constructor
        topo2 = Subzero.TopographyElement(coords)
        @test topo2.coords == coords

        # check when radius is less than  or equal to 0
    end

    @testset "Domain" begin
        # test basic domain
        # domain with wrong directions
        # domain with non-periodic 
        # domain with north < south
        # domain with east < west
    end

    @testset "Floe" begin
        # test generate MC points
        # test with polygon inputs
        # test with coords input
        
    end

    @testset "Model" begin
        # test domain in grid
        # test basic working model
        # test model where domain isn't in grid
        # test model where size of grid and ocean/wind doesn't match
        # test model where types don't match up 
    end
end