using Test, Subzero
import StaticArrays as SA

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
            [CellStresses{Float64}() for i in 1:(g.Nx + 1), j in 1:(g.Ny + 1)],
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

    @testset "Model" begin
        # test domain in grid
        # test basic working model
        # test model where domain isn't in grid
        # test model where size of grid and ocean/atmos doesn't match
        # test model where types don't match up 
    end
end