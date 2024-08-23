@testset "CellStresses" begin
    # Empty CellStresses
    c1 = Subzero.CellStresses()
    @test isempty(c1.τx) && isempty(c1.τy) && isempty(c1.npoints)
    @test typeof(c1) == Subzero.CellStresses{Float64}
    # Non-empty CellStresses
    c2 = Subzero.CellStresses(Float32; npoints = [1, 2], τx = [0.0, 0.0], τy = [4.0, 4.0])
    @test length(c2.npoints) == length(c2.τx) == length(c2.τy) == 2
    @test typeof(c2) == Subzero.CellStresses{Float32}
    # CellStresses with incorrect dimensions
    @test_throws AssertionError Subzero.CellStresses(Float32; npoints = [1], τx = [0.0, 0.0], τy = [4.0, 4.0])
    # Empty CellStresses
    empty!(c2)
    @test isempty(c2.npoints) && isempty(c2.τx) && isempty(c2.τy)
end

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
