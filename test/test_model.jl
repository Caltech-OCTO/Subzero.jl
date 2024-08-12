@testset "Model Creation" begin
    @testset "Ocean" begin
        g = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
        # Large ocean default constructor
        uocn = fill(3.0, g.Nx + 1, g.Ny + 1)
        vocn = fill(4.0, g.Nx + 1, g.Ny + 1)
        tempocn = fill(-2.0, g.Nx + 1, g.Ny + 1)
        τx = fill(0.0, g.Nx + 1, g.Ny + 1)
        τy = τx
        si_frac = fill(0.0, g.Nx + 1, g.Ny + 1)
        hflx_factor = si_frac
        dissolved = si_frac
        ocn = Subzero.Ocean(
            uocn,
            vocn,
            tempocn,
            hflx_factor,
            [IceStressCell{Float64}() for i in 1:(g.Nx + 1), j in 1:(g.Ny + 1)],
            τx,
            τy,
            si_frac,
            dissolved,
            )
        @test ocn.u == uocn
        @test ocn.v == vocn
        @test ocn.temp == tempocn
        @test τx == ocn.τx == ocn.τy
        @test si_frac == ocn.si_frac == ocn.hflx_factor == ocn.dissolved
        @test ocn.u == uocn
        @test ocn.v == vocn
        @test ocn.temp == tempocn
        # Custom constructor
        ocn2 = Subzero.Ocean(g, 3.0, 4.0, -2.0)
        @test ocn.u == ocn2.u
        @test ocn.v == ocn2.v
        @test ocn.temp == ocn2.temp
        @test ocn.si_frac == ocn2.si_frac
        @test ocn.hflx_factor == ocn2.hflx_factor
        @test ocn.τx == ocn2.τx
        @test ocn.τy == ocn2.τy
        # Custom constructor Float32 and Float64
        @test typeof(Subzero.Ocean(Float32, g, 3.0, 4.0, -2.0)) ==
            Subzero.Ocean{Float32}
        @test typeof(Subzero.Ocean(Float64, g, 3.0, 4.0, -2.0)) ==
            Subzero.Ocean{Float64}
    end

    @testset "Atmos" begin
        g = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)

        # Large Atmos default constructor
        uatmos = fill(3.0, g.Nx + 1, g.Ny + 1)
        vatmos = fill(4.0, g.Nx + 1, g.Ny + 1)
        tempatmos = fill(-2.0, g.Nx + 1, g.Ny + 1)
        atmos = Subzero.Atmos(uatmos, vatmos, tempatmos)
        @test atmos.u == uatmos
        @test atmos.v == vatmos
        @test atmos.temp == tempatmos
        # Custom constructor
        atmos2 = Subzero.Atmos(g, 3.0, 4.0, -2.0)
        @test atmos.u == atmos2.u
        @test atmos.v == atmos2.v
        @test atmos.temp == atmos2.temp
        # Custom constructor Float32 anf Float64
        @test typeof(Subzero.Atmos(Float32, g, 3.0, 4.0, -2.0)) ==
            Subzero.Atmos{Float32}
        @test typeof(Subzero.Atmos(Float64, g, 3.0, 4.0, -2.0)) ==
            Subzero.Atmos{Float64}
    end

    @testset "Boundaries" begin
        # Boundaries using BoundaryCoords
        FT = Float64
        g = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
        b1 = Subzero.PeriodicBoundary(North, g)
        b2 = Subzero.OpenBoundary(East, g)
        b3 = Subzero.CollisionBoundary(West, g)
        b4 = Subzero.PeriodicBoundary(South, g)
        b5 = Subzero.MovingBoundary(South, g, 1.0, 2.0)
        @test b1.val == 3e5
        @test typeof(b1) == Subzero.PeriodicBoundary{North, Float64}
        @test b1.coords == [[[-2e5, 3e5], [-2e5, 4.5e5], [6e5, 4.5e5], [6e5, 3e5], [-2e5, 3e5]]]
        @test b2.val == 4e5
        @test typeof(b2) == Subzero.OpenBoundary{East, Float64}
        @test b2.coords == [[[4e5, -1.5e5], [4e5, 4.5e5], [6e5, 4.5e5], [6e5, -1.5e5], [4e5, -1.5e5]]]
        @test b3.val == 0.0
        @test typeof(b3) == Subzero.CollisionBoundary{West, Float64}
        @test b3.coords == [[[-2e5, -1.5e5], [-2e5, 4.5e5], [0.0, 4.5e5], [0.0, -1.5e5], [-2e5, -1.5e5]]]
        @test b4.val == 0.0
        @test typeof(b4) == Subzero.PeriodicBoundary{South, Float64}
        @test b4.coords == [[[-2e5, -1.5e5], [-2e5, 0.0], [6e5, 0.0], [6e5, -1.5e5], [-2e5, -1.5e5]]]
        @test b5.u == 1.0
        @test b5.v == 2.0

        # Test changing CollisionBoundary values and can't change other boundary types
        b5.val = 1.0
        @test b5.val == 1.0
        @test_throws Exception b4.val = 1.0

        # Creation of Float32 and Float64 Boundary
        b32 = Subzero.OpenBoundary(Float32, North, g)
        @test typeof(b32.val) == Float32
        @test typeof(b32) == Subzero.OpenBoundary{North, Float32}
        b64 = Subzero.OpenBoundary(Float64, North, g)
        @test typeof(b64.val) == Float64
        @test typeof(b64) == Subzero.OpenBoundary{North, Float64}

        # Periodic Compat
        @test !Subzero.periodic_compat(b1, b2)
        @test !Subzero.periodic_compat(b2, b1)
        @test Subzero.periodic_compat(b1, b4)
        @test Subzero.periodic_compat(b2, b3)
    end

    @testset "Topography" begin
        coords = [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]]
        poly = Subzero.make_polygon(coords)
        # Polygon Constructor
        topo1 = Subzero.TopographyElement(poly)
        @test topo1.coords == coords
        @test topo1.centroid == [0.5, 0.5]
        @test topo1.rmax == sqrt(0.5)
        topo32 = Subzero.TopographyElement(Float32, poly)
        @test typeof(topo32) == Subzero.TopographyElement{Float32}
        @test typeof(topo32.coords) == Subzero.PolyVec{Float32}
        # Coords Constructor
        topo2 = Subzero.TopographyElement(Float64, coords)
        @test topo2.coords == coords
        @test topo2.centroid == [0.5, 0.5]
        @test topo2.rmax == sqrt(0.5)
        # Basic constructor
        topo3 = TopographyElement(
            Subzero.make_polygon(coords),
            coords,
            [0.5, 0.5],
            sqrt(0.5),
        )
        @test topo3.coords == coords
        # check when radius is less than  or equal to 0
        @test_throws ArgumentError TopographyElement(
            Subzero.make_polygon(coords),
            coords,
            [0.5, 0.5],
            -sqrt(0.5),
        )

        # Create field of topography
        coords_w_hole = [
            [[0.5, 10.0], [0.5, 0.0], [10.0, 0.0], [10.0, 10.0], [0.5, 10.0]],
            [[2.0, 8.0], [2.0, 4.0], [8.0, 4.0], [8.0, 8.0], [2.0, 8.0]]
            ]
        topo_field_64 = initialize_topography_field(
            Float64,
            [coords, coords_w_hole],
        )
        @test length(topo_field_64) == 2
        @test typeof(topo_field_64) <: StructArray{TopographyElement{Float64}}
        @test !Subzero.hashole(topo_field_64.coords[2])

        topo_field_32 = initialize_topography_field(
            Float32,
            [coords, coords_w_hole],
        )
        @test typeof(topo_field_32) <: StructArray{TopographyElement{Float32}}

        @test typeof(initialize_topography_field([coords, coords_w_hole])) <:
            StructArray{TopographyElement{Float64}}
    end

    @testset "Domain" begin
        FT = Float64
        g = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
        b1 = Subzero.PeriodicBoundary(North, g)
        b2 = Subzero.OpenBoundary(East, g)
        b3 = Subzero.CollisionBoundary(West, g)
        b4 = Subzero.PeriodicBoundary(South, g)
        topography = StructArray([Subzero.TopographyElement(
            [[[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]],
        )])
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
        @test_throws ArgumentError Subzero.Domain(
            b1,
            Subzero.OpenBoundary(South, g),
            b2,
            b3,
        )
        @test_throws ArgumentError Subzero.Domain(
            b1,
            b4,
            b2,
            Subzero.PeriodicBoundary(West, g),
        )
        # domain with north < south
        p_placeholder = GI.Polygon([[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.0, 0.0)]])
        @test_throws ArgumentError Subzero.Domain(
            b1,
            Subzero.OpenBoundary(
                South,
                p_placeholder,
                GI.coordinates(p_placeholder),
                6e5,
            ),
            b2,
            b3,
        )
        # domain with east < west
        @test_throws ArgumentError Subzero.Domain(
            b1,
            b4,
            b2,
            Subzero.OpenBoundary(
                West,
                p_placeholder,
                GI.coordinates(p_placeholder),
                6e5,
            ),
        )
    end

    @testset "Model" begin
        # test domain in grid
        # test basic working model
        # test model where domain isn't in grid
        # test model where size of grid and ocean/atmos doesn't match
        # test model where types don't match up 
    end
end