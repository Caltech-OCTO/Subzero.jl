@testset "Model Creation" begin
    # Grid Creation
    @testset "Grid" begin
        # Large grid default constructor
        xlines = collect(-4e5:1e4:4e5)
        ylines = collect(-3e5:1e4:3e5)
        xgrid = repeat(reshape(xlines, 1, :), inner=(length(ylines),1))
        ygrid = repeat(ylines, outer = (1, length(xlines)))
        gbig = Subzero.Grid((61, 81), xlines, ylines, xgrid, ygrid)
        @test gbig.dims == (61, 81)
        @test gbig.xlines == xlines
        @test gbig.ylines == ylines
        @test gbig.xgrid == xgrid
        @test gbig.ygrid == ygrid
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Grid((61, 81), xlines, ylines,
                                                Float64[0.0], ygrid)
        @test_throws ArgumentError Subzero.Grid((61, 81), xlines, ylines,
                                                xgrid, Float64[0.0])
        @test_throws ArgumentError Subzero.Grid((61, 81), Float64[0.0], ylines,
        xgrid, ygrid)
        @test_throws ArgumentError Subzero.Grid((61, 81), xlines, Float64[0.0],
        xgrid, ygrid)
        
        # Non-square grid using custom constructor
        g1 = Subzero.Grid(10, 8, 2)
        @test g1.dims == (9, 11)
        @test g1.xlines == collect(-10.0:2:10.0)
        @test g1.ylines == collect(-8:2:8)
        @test g1.xgrid[1, :] == g1.xlines
        @test g1.ygrid[:, 1] == g1.ylines
        @test typeof(g1) == Subzero.Grid{Float64}
        # Uneven grid size creation (grid cut short) using custom constructor
        g2 = Subzero.Grid(10.5, 8.0, 2.5)
        @test g2.dims == (7, 9)
        @test g2.xlines == collect(-10.5:2.5:10.5)
        @test size(g2.xgrid) == g2.dims
        # Custom constructor Float32
        @test typeof(Subzero.Grid(10, 8, 2, Float32)) == Subzero.Grid{Float32}
    end

    @testset "Ocean" begin
        g = Subzero.Grid(4e5, 3e5, 1e4)
        # Large ocean default constructor
        uocn = fill(3.0, g.dims)
        vocn = fill(4.0, g.dims)
        tempocn = fill(-2.0, g.dims)
        ocn = Subzero.Ocean(uocn, vocn, tempocn)
        @test ocn.u == uocn
        @test ocn.v == vocn
        @test ocn.temp == tempocn
        # Default constructor fails for non-matching dimensions
        @test_throws ArgumentError Subzero.Ocean(Float64[0.0], vocn, tempocn)
        @test_throws ArgumentError Subzero.Ocean(uocn, Float64[0.0], tempocn)
        @test_throws ArgumentError Subzero.Ocean(uocn, vocn, Float64[0.0])
        # Custom constructor
        ocn2 = Subzero.Ocean(g, 3.0, 4.0, -2.0)
        @test ocn.u == ocn2.u
        @test ocn.v == ocn2.v
        @test ocn.temp == ocn2.temp
        # Custom constructor Float32
        @test typeof(Subzero.Ocean(g, 3.0, 4.0, -2.0, Float32)) ==
              Subzero.Ocean{Float32}
    end

    @testset "Wind" begin
        g = Subzero.Grid(4e5, 3e5, 1e4)
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
        # Periodic Compat
        @test !Subzero.periodic_compat(b1, b2)
        @test !Subzero.periodic_compat(b2, b1)
        @test Subzero.periodic_compat(b1, b4)
        @test Subzero.periodic_compat(b2, b3)
    end


    # Ocean Creation

    # Wind Creation

    # Boundary Creation

    # Domain Creation

    # Topography Creation

    # Floe Creation

    # Model Creation
end