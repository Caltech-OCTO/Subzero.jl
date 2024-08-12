using Test, Subzero

# # Default constructor fails for non-matching dimensions
# @test_throws ArgumentError Subzero.RegRectilinearGrid(
#     80,
#     50,
#     0,
#     1e5,
#     0,
#     1e5,
#     1e3,  # 1e5/80 = 1.25e3 ≂̸ 1e3
#     2e3,
#     [CellFloes{Float64}() for i in 1:81, j in 1:51],
# )
# @test_throws ArgumentError Subzero.RegRectilinearGrid(
#     80,
#     50,
#     0,
#     1e5,
#     0,
#     1e5,
#     1.25e3,
#     1e3,  # 1e5/50 = 2e3 ≂̸ 1e3
#     [CellFloes{Float64}() for i in 1:81, j in 1:51],
# )
# @test_throws ArgumentError Subzero.RegRectilinearGrid(
#     80,
#     50,
#     0,
#     1e5,
#     0,
#     1e5,
#     1.25e3,
#     2e3,
#     [CellFloes{Float64}() for i in 1:71, j in 1:51],  # wrong dims
# )

# Non-square grid using constructor with Δx and Δy
g1 = Subzero.RegRectilinearGrid(; x0 = -10, xf = 10, y0 = -8, yf = 8, Δx = 2, Δy = 4)
@test g1.Nx == 10
@test g1.Ny == 4
@test g1.x0 == -10
@test g1.xf == 10
@test g1.y0 == -8
@test g1.yf == 8
@test size(g1.floe_locations) == (11, 5)
@test typeof(g1) == Subzero.RegRectilinearGrid{Float64}

# Non-square grid using constructor with Nx and Ny
g2 = Subzero.RegRectilinearGrid(; x0 = -10, xf = 10, y0 = -8, yf = 8, Nx = 10, Ny = 4)
@test g2.x0 == -10
@test g2.xf == 10
@test g2.y0 == -8
@test g2.yf == 8
@test g2.Δx == 2
@test g2.Δy == 4
@test size(g2.floe_locations) == (11, 5)
@test typeof(g2) == Subzero.RegRectilinearGrid{Float64}

# Custom constructor Float32 and Float64
@test typeof(Subzero.RegRectilinearGrid(Float32; x0 = 0, xf = 10, y0 = 0, yf = 8, Δx = 2, Δy = 2)) == Subzero.RegRectilinearGrid{Float32}
@test typeof(Subzero.RegRectilinearGrid(Float64; x0 = 0, xf = 10, y0 = 0, yf = 8, Δx = 2, Δy = 2)) == Subzero.RegRectilinearGrid{Float64}
