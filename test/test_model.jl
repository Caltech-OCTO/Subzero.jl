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
        si_frac = τx
        ocn = Subzero.Ocean(uocn, vocn, tempocn, τx, τy, si_frac)
        @test ocn.u == uocn
        @test ocn.v == vocn
        @test ocn.temp == tempocn
        @test ocn.τx == τx == ocn.τy == ocn.si_frac
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Ocean(Float64[0.0], vocn, tempocn, τx, τy, si_frac)
        @test_throws ArgumentError Subzero.Ocean(uocn, Float64[0.0], tempocn, τx, τy, si_frac)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, Float64[0.0], τx, τy, si_frac)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, Float64[0.0], τy, si_frac)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, tempocn, τx, τy, Float64[0.0])
        # Custom constructor
        ocn2 = Subzero.Ocean(g, 3.0, 4.0, -2.0)
        @test ocn.u == ocn2.u
        @test ocn.v == ocn2.v
        @test ocn.temp == ocn2.temp
        @test ocn.τx == ocn2.τx
        @test ocn.si_frac == ocn2.si_frac
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
        b1 = Subzero.Boundary(Subzero.PeriodicBC(), 0.0)
        b2 = Subzero.Boundary(Subzero.OpenBC(), 1.0)
        b3 = Subzero.Boundary(Subzero.CollisionBC(), 2.0)
        b4 = Subzero.Boundary(Subzero.PeriodicBC(), 5.0)
        @test b1.bc == Subzero.PeriodicBC()
        @test b1.val == 0.0
        b32 = Subzero.Boundary(Subzero.OpenBC(), 1, Float32)
        @test typeof(b32.val) == Float32

        # Periodic Compat
        @test !Subzero.periodic_compat(b1, b2)
        @test !Subzero.periodic_compat(b2, b1)
        @test Subzero.periodic_compat(b1, b4)
        @test Subzero.periodic_compat(b2, b3)

        # RectangularDomain default constructor
        rdomain1 = Subzero.RectangleDomain(b4, b1, b3, b2)
        @test rdomain1.north == b4
        @test rdomain1.south == b1
        @test rdomain1.east == b3

        # RectangularDomain fails
        @test_throws ArgumentError Subzero.RectangleDomain(b1, b4, b3, b2)
        @test_throws ArgumentError Subzero.RectangleDomain(b4, b1, b2, b3)
        @test_throws ArgumentError Subzero.RectangleDomain(b4, b1, b2, b2)
        @test_throws ArgumentError Subzero.RectangleDomain(b4, b2, b2, b3)
        @test_throws MethodError Subzero.RectangleDomain(b4, b1, b3, b32)

        # RectangularDomain Grid constructor
        g = Subzero.Grid(4e5, 3e5, 1e4, 1e4)
        g32 = Subzero.Grid(4e5, 3e5, 1e4, 1e4, Float32)
        rdomain2 = Subzero.RectangleDomain(g)
        @test rdomain2.north == Subzero.Boundary(Subzero.OpenBC(), 3e5)
        @test rdomain2.south == Subzero.Boundary(Subzero.OpenBC(), 0.0)
        rdomain3 = Subzero.RectangleDomain(g, northBC = Subzero.CollisionBC(),
                                           eastBC = Subzero.CompressionBC())
        @test rdomain3.north.bc == Subzero.CollisionBC()
        @test rdomain3.east.bc == Subzero.CompressionBC()
        @test typeof(Subzero.RectangleDomain(g32).north.val) == Float32

        # CircleDomain default constructor
        cdomain1 = Subzero.CircleDomain(5.0, [1.0, 1.0], Subzero.OpenBC())
        @test cdomain1.radius == 5.0
        @test cdomain1.centroid == [1.0, 1.0]
        @test cdomain1.bc == Subzero.OpenBC()

        # CircleDomain fails
        @test_throws ArgumentError Subzero.CircleDomain(-5.0, [1.0, 1.0],
                                                        Subzero.OpenBC())
        # CircleDomain Grid constructor
        cdomain2 = Subzero.CircleDomain(g)
        @test cdomain2.radius == 3e5/2
        @test cdomain2.centroid == [2e5, 1.5e5]
        @test cdomain2.bc == Subzero.OpenBC()
        cdomain3 = Subzero.CircleDomain(g, bc = Subzero.CollisionBC())
        @test cdomain3.bc == Subzero.CollisionBC()

        # Grid from domains
        rdomain_grid = Subzero.Grid(rdomain2, (30, 40))
        @test rdomain_grid.dims == g.dims
        @test rdomain_grid.xg == g.xg
        @test rdomain_grid.xc == g.xc
        cdomain_grid = Subzero.Grid(cdomain2, (30, 40))
        @test cdomain_grid.dims == g.dims
        @test cdomain_grid.xg == collect(5e4:7.5e3:3.5e5)
        @test cdomain_grid.xc == collect(5.375e4:7.5e3:3.4625e5)
    end

    @testset "Topography" begin
        coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]]
        poly = LibGEOS.Polygon(coords)
        cent = [0.5, 0.5]
        # Polygon Constructor
        topo1 = Subzero.Topography(poly, 0.5)
        @test topo1.centroid == cent
        @test topo1.coords == coords
        @test topo1.height == 0.5
        @test topo1.area == 1.0
        @test topo1.rmax == sqrt(0.5^2 + 0.5^2)
        topo32 = Subzero.Topography(poly, 0.5, Float32)
        @test typeof(topo32) == Subzero.Topography{Float32}
        @test typeof(topo32.coords) == Subzero.PolyVec{Float32}
        # Coords Constructor
        topo2 = Subzero.Topography(coords, 0.5)
        @test topo2.height == 0.5
        @test topo2.area == 1.0
        # Constructor fails
        @test_throws ArgumentError Subzero.Topography(poly, -1.0)
        @test_throws ArgumentError Subzero.Topography(cent, coords, 1., -1., .5)
        @test_throws ArgumentError Subzero.Topography(cent, coords, 1., 1., -.5)

    end

    @testset "Floe" begin
        
    end

    @testset "Model" begin
        
    end
end