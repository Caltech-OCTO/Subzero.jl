@testset "Model Creation" begin
    # Grid Creation
    @testset "Grid" begin
        # Large grid default constructor
        xg = collect(-4e5:1e4:4e5)
        yg = collect(-3e5:1e4:3e5)
        xc = collect(-3.95e5:1e4:3.95e5)
        yc = collect(-2.95e5:1e4:2.95e5)
        gbig = Subzero.Grid((length(yc), length(xc)), xg, yg, xc, yc)
        @test gbig.dims == (60, 80)
        @test gbig.xg == xg
        @test gbig.yg == yg
        @test gbig.xc == xc
        @test gbig.yc == yc
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Grid((60, 80), xg, yg,
                                                Float64[0.0], yc)
        @test_throws ArgumentError Subzero.Grid((60, 80), xg, yg,
                                                xc, Float64[0.0])
        @test_throws ArgumentError Subzero.Grid((60, 80), Float64[0.0], yg,
                                                xc, yc)
        @test_throws ArgumentError Subzero.Grid((60, 80), xg, Float64[0.0],
                                                xc, yc)
        
        # Non-square grid using custom constructor
        g1 = Subzero.Grid(-10, 10, -8, 8, 2, 4)
        @test g1.dims == (4, 10)
        @test g1.xg == collect(-10.0:2:10.0)
        @test g1.yg == collect(-8:4:8)
        @test g1.xc == collect(-9.0:2:9.0)
        @test g1.yc == collect(-6.0:4:6.0)
        @test typeof(g1) == Subzero.Grid{Float64}
        # Uneven grid size creation (grid cut short) using custom constructor
        g2 = Subzero.Grid(10.5, 8.0, 2.5, 2)
        @test g2.dims == (4, 4)
        @test g2.xg == collect(0.0:2.5:10.5)
        @test g2.yg == collect(0.0:2.0:8.0)
        @test g2.xc == collect(1.25:2.5:8.75)
        @test g2.yc == collect(1.0:2.0:7.0)
        # Custom constructor Float32
        @test typeof(Subzero.Grid(10, 8, 2, 2, Float32)) ==
              Subzero.Grid{Float32}
        # Grid constructor with dimensions
        g3 = Subzero.Grid(-10, 10, -8, 8, (4, 10))
        @test g3.dims == (4,10)
        @test g3.xg == g1.xg
        @test g3.yc == g1.yc
    end

    @testset "Ocean" begin
        g = Subzero.Grid(4e5, 3e5, 1e4, 1e4)
        # Large ocean default constructor
        uocn = fill(3.0, g.dims)
        vocn = fill(4.0, g.dims)
        tempocn = fill(-2.0, g.dims)
        τx = fill(0.0, g.dims)
        τy = τx
        si_area = τx
        hflx = τx
        ocn = Subzero.Ocean(uocn, vocn, tempocn, hflx, τx, τy, si_area)
        @test ocn.u == uocn
        @test ocn.v == vocn
        @test ocn.temp == tempocn
        @test ocn.τx == τx == ocn.τy == ocn.si_area == ocn.hflx
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Ocean(Float64[0.0], vocn, tempocn, hflx, τx, τy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, Float64[0.0], tempocn, hflx, τx, τy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, Float64[0.0], hflx, τx, τy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, Float64[0.0], τx, τy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, hflx, Float64[0.0], τy, si_area)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, hflx, τx, τy, Float64[0.0])
        # Custom constructor
        ocn2 = Subzero.Ocean(g, 3.0, 4.0, -2.0)
        @test ocn.u == ocn2.u
        @test ocn.v == ocn2.v
        @test ocn.temp == ocn2.temp
        @test ocn.τx == ocn2.τx
        @test ocn.si_frac == ocn2.si_frac
        @test ocn.hflx == ocn2.hflx
        # Custom constructor Float32
        @test typeof(Subzero.Ocean(g, 3.0, 4.0, -2.0, Float32)) ==
              Subzero.Ocean{Float32}
    end

    @testset "Wind" begin
        g = Subzero.Grid(4e5, 3e5, 1e4, 1e4)
        # Large wind default constructor
        uwind = fill(3.0, g.dims)
        vwind = fill(4.0, g.dims)
        tempwind = fill(-2.0, g.dims)
        wind = Subzero.Wind(uwind, vwind, tempwind)
        @test wind.u == uwind
        @test wind.v == vwind
        @test wind.temp == tempwind
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Wind(Float64[0.0], vwind, tempwind)
        @test_throws ArgumentError Subzero.Wind(uwind, Float64[0.0], tempwind)
        @test_throws ArgumentError Subzero.Wind(uwind, vwind, Float64[0.0])
        # Custom constructor
        wind2 = Subzero.Wind(g, 3.0, 4.0, -2.0)
        @test wind.u == wind2.u
        @test wind.v == wind2.v
        @test wind.temp == wind2.temp
        # Custom constructor Float32
        @test typeof(Subzero.Wind(g, 3.0, 4.0, -2.0, Float32)) ==
              Subzero.Wind{Float32}
    end

    @testset "Boundary and Domain" begin
        # Boundary
        g = Subzero.Grid(4e5, 3e5, 1e4, 1e4)
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

        # Need to add in ability to run with Float32
        #b32 = Subzero.OpenBoundary(Subzero.OpenBC(), 1, Float32)
        #@test typeof(b32.val) == Float32

        # Periodic Compat
        @test !Subzero.periodic_compat(b1, b2)
        @test !Subzero.periodic_compat(b2, b1)
        @test Subzero.periodic_compat(b1, b4)
        @test Subzero.periodic_compat(b2, b3)

        # RectangularDomain default constructor
        rdomain1 = Subzero.Domain(b1, b4, b2, b3)
        @test rdomain1.north == b1

        # RectangularDomain fails
        # need checks for val not being correct
        @test_throws ArgumentError Subzero.Domain(b4, b2, b2, b3)
        #@test_throws MethodError Subzero.Domain(b4, b1, b3, b32)
    end

    @testset "Topography" begin
        coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]]
        poly = LibGEOS.Polygon(coords)
        # Polygon Constructor
        topo1 = Subzero.Topography(poly)
        @test topo1.coords == coords
        topo32 = Subzero.Topography(poly, Float32)
        @test typeof(topo32) == Subzero.Topography{Float32}
        @test typeof(topo32.coords) == Subzero.PolyVec{Float32}
        # Coords Constructor
        topo2 = Subzero.Topography(coords)
        @test topo2.coords == coords
    end

    @testset "Floe" begin
        
    end

    @testset "Model" begin
        
    end
end