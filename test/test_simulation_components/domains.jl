using Test, Subzero
import GeometryOps.GeoInterface as GI

    @testset "Directions" begin
        FT = Float64
        x0, xf, y0, yf = 0.0, 1e5, -5e4, 5e4
        Δx, Δy = (xf - x0) / 2, (yf - y0) / 2 
        grid = RegRectilinearGrid(; x0, xf, y0, yf, Nx = 10, Ny = 10)
        # North boundary
        n_poly, n_val = Subzero._grid_boundary_info(FT, North, grid)
        n_point_set = Set(GI.getpoint(n_poly))
        @test issetequal(n_point_set, Set(((x0 - Δx, yf), (x0 - Δx, yf + Δy), (xf + Δx, yf + Δy), (xf + Δx, yf))))
        @test n_val == grid.yf
        # South boundary
        s_poly, s_val = Subzero._grid_boundary_info(FT, South, grid)
        s_point_set = Set(GI.getpoint(s_poly))
        @test issetequal(s_point_set, Set(((x0 - Δx, y0), (x0 - Δx, y0 - Δy), (xf + Δx, y0 - Δy), (xf + Δx, y0))))
        @test s_val == grid.y0
        # East boundary
        e_poly, e_val = Subzero._grid_boundary_info(FT, East, grid)
        e_point_set = Set(GI.getpoint(e_poly))
        @test issetequal(e_point_set, Set(((xf, y0 - Δy), (xf, yf + Δy), (xf + Δx, y0 - Δy), (xf + Δx, yf + Δy))))
        @test e_val == grid.xf
        # West boundary
        w_poly, w_val = Subzero._grid_boundary_info(FT, West, grid)
        w_point_set = Set(GI.getpoint(w_poly))
        @test issetequal(w_point_set, Set(((x0, y0 - Δy), (x0, yf + Δy), (x0 - Δx, y0 - Δy), (x0 - Δx, yf + Δy))))
        @test w_val == grid.x0
    end

    # Boundaries using BoundaryCoords
    FT = Float64
    g = Subzero.RegRectilinearGrid(; x0 = 0, xf = 4e5, y0 = 0, yf = 3e5, Δx = 1e4, Δy = 1e4)
    b1 = Subzero.PeriodicBoundary(North; grid = g)
    b2 = Subzero.OpenBoundary(East; grid = g)
    b3 = Subzero.CollisionBoundary(West; grid = g)
    b4 = Subzero.PeriodicBoundary(South; grid = g)
    b5 = Subzero.MovingBoundary(South; grid = g, u = 1.0, v = 2.0)
    @test b1.val == 3e5
    @test typeof(b1) == Subzero.PeriodicBoundary{North, Float64}
    # @test b1.coords == [[[-2e5, 3e5], [-2e5, 4.5e5], [6e5, 4.5e5], [6e5, 3e5], [-2e5, 3e5]]]
    @test b2.val == 4e5
    @test typeof(b2) == Subzero.OpenBoundary{East, Float64}
    # @test b2.coords == [[[4e5, -1.5e5], [4e5, 4.5e5], [6e5, 4.5e5], [6e5, -1.5e5], [4e5, -1.5e5]]]
    @test b3.val == 0.0
    @test typeof(b3) == Subzero.CollisionBoundary{West, Float64}
    # @test b3.coords == [[[-2e5, -1.5e5], [-2e5, 4.5e5], [0.0, 4.5e5], [0.0, -1.5e5], [-2e5, -1.5e5]]]
    @test b4.val == 0.0
    @test typeof(b4) == Subzero.PeriodicBoundary{South, Float64}
    # @test b4.coords == [[[-2e5, -1.5e5], [-2e5, 0.0], [6e5, 0.0], [6e5, -1.5e5], [-2e5, -1.5e5]]]
    @test b5.u == 1.0
    @test b5.v == 2.0

    # Test changing CollisionBoundary values and can't change other boundary types
    b5.val = 1.0
    @test b5.val == 1.0
    @test_throws Exception b4.val = 1.0

    # Creation of Float32 and Float64 Boundary
    b32 = Subzero.OpenBoundary(North, Float32; grid = g)
    @test typeof(b32.val) == Float32
    @test typeof(b32) == Subzero.OpenBoundary{North, Float32}
    b64 = Subzero.OpenBoundary(North, Float64; grid = g)
    @test typeof(b64.val) == Float64
    @test typeof(b64) == Subzero.OpenBoundary{North, Float64}

    # Periodic Compat
    @test !Subzero.periodic_compat(b1, b2)
    @test !Subzero.periodic_compat(b2, b1)
    @test Subzero.periodic_compat(b1, b4)
    @test Subzero.periodic_compat(b2, b3)


    # Test compression boundaries movement - TODO: Move these to different file
    # Lx = 1e5
    # Ly = Lx
    # grid = RegRectilinearGrid(; x0 = -Lx, xf = Lx, y0 = -Ly, yf = Ly, Δx = 1e4, Δy = 1e4)

    # nc_boundary = MovingBoundary(North; grid, 0.0, -0.1)
    # nc_coords = deepcopy(nc_boundary.coords)
    # sc_boundary = MovingBoundary(South; grid, 0.0, 0.1)
    # sc_coords = deepcopy(sc_boundary.coords)
    # ec_boundary = MovingBoundary(East; grid, 0.1, 0.0)
    # ec_coords = deepcopy(ec_boundary.coords)
    # wc_boundary = MovingBoundary(West; grid, 0.1, 0.0)
    # wc_coords = deepcopy(wc_boundary.coords)
    # cdomain = Domain(nc_boundary, sc_boundary, ec_boundary, wc_boundary)
    # Subzero.update_boundaries!(cdomain, 10)
    # Subzero.translate!(nc_coords, 0, -1)
    # @test nc_coords == nc_boundary.coords
    # @test nc_boundary.val == 1e5 - 1
    # Subzero.translate!(sc_coords, 0, 1)
    # @test sc_coords == sc_boundary.coords
    # @test sc_boundary.val == -1e5 + 1
    # Subzero.translate!(ec_coords, 1, 0)
    # @test ec_coords == ec_boundary.coords
    # @test ec_boundary.val == 1e5 + 1
    # Subzero.translate!(wc_coords, 1, 0)
    # @test wc_coords == wc_boundary.coords
    # @test wc_boundary.val == -1e5 + 1


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