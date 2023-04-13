@testset "Model Creation" begin
    # Grid Creation
    @testset "Grid" begin
        # Large grid default constructor
        xg = collect(-4e5:1e4:4e5)
        yg = collect(-3e5:1e4:3e5)
        xc = collect(-3.95e5:1e4:3.95e5)
        yc = collect(-2.95e5:1e4:2.95e5)
        gbig = Subzero.RegRectilinearGrid(
            (length(yc), length(xc)),
            xg,
            yg,
            xc,
            yc,
            [CellFloes{Float64}() for i in 1:61, j in 1:81]
        )
        @test gbig.dims == (60, 80)
        @test gbig.xg == xg
        @test gbig.yg == yg
        @test gbig.xc == xc
        @test gbig.yc == yc
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.RegRectilinearGrid(
            (60, 80),
            xg,
            yg,
            Float64[0.0],
            yc,
            [CellFloes{Float64}() for i in 1:61, j in 1:81],
        )
        @test_throws ArgumentError Subzero.RegRectilinearGrid(
            (60, 80),
            xg,
            yg,
            xc,
            Float64[0.0],
            [CellFloes{Float64}() for i in 1:61, j in 1:81],
        )
        @test_throws ArgumentError Subzero.RegRectilinearGrid(
            (60, 80),
            Float64[0.0],
            yg,
            xc,
            yc,
            [CellFloes{Float64}() for i in 1:61, j in 1:81],
        )
        @test_throws ArgumentError Subzero.RegRectilinearGrid(
            (60, 80),
            xg,
            Float64[0.0],
            xc,
            yc,
            [CellFloes{Float64}() for i in 1:61, j in 1:81],
        )
        
        # Non-square grid using custom constructor
        g1 = Subzero.RegRectilinearGrid(
            Float64,
            (-10, 10),
            (-8, 8),
            2,
            4,
        )
        @test g1.dims == (4, 10)
        @test g1.xg == collect(-10.0:2:10.0)
        @test g1.yg == collect(-8:4:8)
        @test g1.xc == collect(-9.0:2:9.0)
        @test g1.yc == collect(-6.0:4:6.0)
        @test typeof(g1) == Subzero.RegRectilinearGrid{Float64}
        # Uneven grid size creation (grid cut short) using custom constructor
        g2 = Subzero.RegRectilinearGrid(
            Float64,
            (0.0, 10.5),
            (0.0, 8.0),
            2.5,
            2,
        )
        @test g2.dims == (4, 4)
        @test g2.xg == collect(0.0:2.5:10.5)
        @test g2.yg == collect(0.0:2.0:8.0)
        @test g2.xc == collect(1.25:2.5:8.75)
        @test g2.yc == collect(1.0:2.0:7.0)
        # Custom constructor Float32
        @test typeof(Subzero.RegRectilinearGrid(
            Float32,
            (0, 10),
            (0, 8),
            2,
            2,
        )) == Subzero.RegRectilinearGrid{Float32}
        # Grid constructor with dimensions
        g3 = Subzero.RegRectilinearGrid(
            Float64,
            (-10, 10),
            (-8, 8),
            (4, 10),
        )
        @test g3.dims == (4,10)
        @test g3.xg == g1.xg
        @test g3.yc == g1.yc
    end

    @testset "Ocean" begin
        g = Subzero.RegRectilinearGrid(
            Float64,
            (0, 4e5),
            (0, 3e5),
            1e4,
            1e4,
        )
        # Large ocean default constructor
        uocn = fill(3.0, g.dims .+ 1)
        vocn = fill(4.0, g.dims .+ 1)
        tempocn = fill(-2.0, g.dims .+ 1)
        fx = fill(0.0, g.dims .+ 1)
        fy = fx
        τx = fx
        τy = fx
        si_frac = fx
        hflx_factor = fx
        ocn = Subzero.Ocean(
            uocn,
            vocn,
            tempocn,
            hflx_factor,
            [IceStressCell{Float64}() for i in 1:g.dims[1] + 1, j in 1:g.dims[2] + 1],
            τx,
            τy,
            si_frac,
        )
        @test ocn.u == uocn
        @test ocn.v == vocn
        @test ocn.temp == tempocn
        @test ocn.si_frac == ocn.hflx_factor == ocn.τx == ocn.τx
        # Custom constructor
        ocn2 = Subzero.Ocean(Float64, g, 3.0, 4.0, -2.0)
        @test ocn.u == ocn2.u
        @test ocn.v == ocn2.v
        @test ocn.temp == ocn2.temp
        @test ocn.si_frac == ocn2.si_frac
        @test ocn.hflx_factor == ocn2.hflx_factor
        @test ocn.τx == ocn2.τx
        @test ocn.τy == ocn2.τy
        # Custom constructor Float32
        @test typeof(Subzero.Ocean(Float32, g, 3.0, 4.0, -2.0)) ==
              Subzero.Ocean{Float32}
    end

    @testset "Atmos" begin
        g = Subzero.RegRectilinearGrid(
            Float64,
            (0, 4e5),
            (0, 3e5),
            1e4,
            1e4,
        )
        # Large Atmos default constructor
        uatmos = fill(3.0, g.dims .+ 1)
        vatmos = fill(4.0, g.dims .+ 1)
        tempatmos = fill(-2.0, g.dims .+ 1)
        atmos = Subzero.Atmos(uatmos, vatmos, tempatmos)
        @test atmos.u == uatmos
        @test atmos.v == vatmos
        @test atmos.temp == tempatmos
        # Custom constructor
        atmos2 = Subzero.Atmos(g, 3.0, 4.0, -2.0)
        @test atmos.u == atmos2.u
        @test atmos.v == atmos2.v
        @test atmos.temp == atmos2.temp
        # Custom constructor Float32
        @test typeof(Subzero.Atmos(g, 3.0, 4.0, -2.0, Float32)) ==
              Subzero.Atmos{Float32}
    end

    @testset "Boundaries" begin
        # Boundaries using BoundaryCoords
        g = Subzero.RegRectilinearGrid(
            Float64,
            (0, 4e5),
            (0, 3e5),
            1e4,
            1e4,
        )
        b1 = Subzero.PeriodicBoundary(g, Subzero.North())
        b2 = Subzero.OpenBoundary(g, Subzero.East())
        b3 = Subzero.CollisionBoundary(g, Subzero.West())
        b4 = Subzero.PeriodicBoundary(g, Subzero.South())
        b5 = Subzero.CompressionBoundary(g, Subzero.South(), 1.0)
        @test b1.val == 3e5
        @test typeof(b1) == Subzero.PeriodicBoundary{Subzero.North, Float64}
        @test b1.coords == [[[-2e5, 3e5], [-2e5, 4.5e5], [6e5, 4.5e5], [6e5, 3e5], [-2e5, 3e5]]]
        @test b2.val == 4e5
        @test typeof(b2) == Subzero.OpenBoundary{Subzero.East, Float64}
        @test b2.coords == [[[4e5, -1.5e5], [4e5, 4.5e5], [6e5, 4.5e5], [6e5, -1.5e5], [4e5, -1.5e5]]]
        @test b3.val == 0.0
        @test typeof(b3) == Subzero.CollisionBoundary{Subzero.West, Float64}
        @test b3.coords == [[[-2e5, -1.5e5], [-2e5, 4.5e5], [0.0, 4.5e5], [0.0, -1.5e5], [-2e5, -1.5e5]]]
        @test b4.val == 0.0
        @test typeof(b4) == Subzero.PeriodicBoundary{Subzero.South, Float64}
        @test b4.coords == [[[-2e5, -1.5e5], [-2e5, 0.0], [6e5, 0.0], [6e5, -1.5e5], [-2e5, -1.5e5]]]
        @test b5.velocity == 1.0
        @test b5.val == 0.0

        # Test changing CollisionBoundary values and can't change other boundary types
        b5.val = 1.0
        @test b5.val == 1.0
        @test_throws Exception b4.val = 1.0

        # Creation of Float32 Boundary
        b32 = Subzero.OpenBoundary(g, Subzero.North(), Float32)
        @test typeof(b32.val) == Float32
        @test typeof(b32) == Subzero.OpenBoundary{Subzero.North, Float32}

        # Periodic Compat
        @test !Subzero.periodic_compat(b1, b2)
        @test !Subzero.periodic_compat(b2, b1)
        @test Subzero.periodic_compat(b1, b4)
        @test Subzero.periodic_compat(b2, b3)
    end

    @testset "Topography" begin
        coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]]
        poly = LG.Polygon(coords)
        # Polygon Constructor
        topo1 = Subzero.TopographyElement(poly)
        @test topo1.coords == coords
        @test topo1.centroid == [0.5, 0.5]
        @test topo1.rmax == sqrt(0.5)
        topo32 = Subzero.TopographyElement(poly, Float32)
        @test typeof(topo32) == Subzero.TopographyElement{Float32}
        @test typeof(topo32.coords) == Subzero.PolyVec{Float32}
        # Coords Constructor
        topo2 = Subzero.TopographyElement(coords)
        @test topo2.coords == coords
        @test topo2.centroid == [0.5, 0.5]
        @test topo2.rmax == sqrt(0.5)
        # Basic constructor
        topo3 = TopographyElement(coords, [0.5, 0.5], sqrt(0.5))
        @test topo3.coords == coords
        # check when radius is less than  or equal to 0
        @test_throws ArgumentError TopographyElement(coords, [0.5, 0.5], -sqrt(0.5))
    end

    @testset "Domain" begin
        g = Subzero.RegRectilinearGrid(
            Float64,
            (0, 4e5),
            (0, 3e5),
            1e4,
            1e4,
        )
        b1 = Subzero.PeriodicBoundary(g, Subzero.North())
        b2 = Subzero.OpenBoundary(g, Subzero.East())
        b3 = Subzero.CollisionBoundary(g, Subzero.West())
        b4 = Subzero.PeriodicBoundary(g, Subzero.South())
        topography = StructArray([Subzero.TopographyElement([[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]])])
        # test basic domain with no topography
        rdomain1 = Subzero.Domain(b1, b4, b2, b3)
        @test rdomain1.north == b1
        @test rdomain1.south == b4
        @test rdomain1.east == b2
        @test rdomain1.west == b3
        @test isempty(rdomain1.topography)
        # test basic domain with topography
        rdomain2 = Subzero.Domain(b1, b4, b2, b3, topography)
        @test isempty(rdomain1.topography)
        # domain with wrong directions
        @test_throws MethodError Subzero.Domain(b4, b2, b2, b3)
        # domain with non-periodic 
        @test_throws ArgumentError Subzero.Domain(b1, Subzero.OpenBoundary(g, Subzero.South()), b2, b3)
        @test_throws ArgumentError Subzero.Domain(b1, b4, b2, Subzero.PeriodicBoundary(g, Subzero.West()))
        # domain with north < south
        @test_throws ArgumentError Subzero.Domain(b1, Subzero.OpenBoundary(Subzero.PolyVec{Float64}(undef, 0), 6e5, Subzero.South()), b2, b3)
        # domain with east < west
        @test_throws ArgumentError Subzero.Domain(b1, b4, b2, Subzero.OpenBoundary(Subzero.PolyVec{Float64}(undef, 0), 6e5, Subzero.West()))
    end

    @testset "Model" begin
        # test domain in grid
        # test basic working model
        # test model where domain isn't in grid
        # test model where size of grid and ocean/atmos doesn't match
        # test model where types don't match up 
    end
end