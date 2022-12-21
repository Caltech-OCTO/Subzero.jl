using Subzero, StructArrays, Statistics, JLD2
import LibGEOS as LG

# User Inputs
const T = Float64::DataType

const Lx = 7.5e4
const Ly = 7.5e4
const Δgrid = 10000
const h_mean = 0.25
const Δh = 0.0
const Δt = 20
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

# Model instantiation
grid = RegRectilinearGrid(-Lx, Lx, -Ly, Ly, Δgrid, Δgrid)

uocn = vcat(hcat(0.5ones(8,8), zeros(8, 8)), hcat(zeros(8,8), -ones(8,8)))
vocn = vcat(hcat(zeros(8,8), 0.5ones(8, 8)), hcat(zeros(8,8), -ones(8,8)))
temp = zeros(16, 16)
ocean = Ocean(uocn, vocn, temp)

wind = Wind(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), fill(-20.0, grid.dims .+ 1))

# Domain creation - boundaries and topography
nboundary = OpenBoundary(grid, North())
sboundary = OpenBoundary(grid, South())
eboundary = OpenBoundary(grid, East())
wboundary = OpenBoundary(grid, West())

domain = Subzero.Domain(nboundary, sboundary, eboundary, wboundary)
\
# Floe instantiation
nfloes = 100
file = jldopen("examples/floe_shapes.jld2", "r")
nfloes = nfloes > size(file["floe_vertices"], 1) ? size(file["floe_vertices"], 1) : nfloes
floe_arr = vec(StructArray([Floe(vert, h_mean, Δh) for vert in file["floe_vertices"][1:nfloes]]))
close(file)

model = Model(grid, ocean, wind, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
#consts = Constants(E = modulus, Cd_io = 0.0, Cd_ia = 0.0, Cd_ao = 0.0, f = 0.0, μ = 0.0)  # collisions without friction 
simulation = Simulation(model = model, consts = consts, Δt = Δt, nΔt = 3000, COLLISION = true)

# Output setup
gridwriter = GridOutputWriter([GridOutput(i) for i in 1:9], 10, "g.nc", grid, (10, 10))
floewriter = FloeOutputWriter([FloeOutput(i) for i in [3:4; 6; 9:11; 22:26]], 50, "f.nc", grid)

# Run simulation
run!(simulation, [floewriter])