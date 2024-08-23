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
