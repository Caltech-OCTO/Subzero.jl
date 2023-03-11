using Subzero, StructArrays, Statistics, JLD2
import LibGEOS as LG

# User Inputs
const FT = Float64

const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const Δt = 20
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

# Model instantiation
grid = RegRectilinearGrid(
    FT,
    (-Lx, Lx),
    (-Ly, Ly),
    Δgrid,
    Δgrid,
)

ocean = Ocean(grid, -0.2, 0.0, -1.0)

atmos = Atmos(grid, 0.0, 0.0, -3.0)

# Domain creation - boundaries and topography
nboundary = OpenBoundary(grid, North())
sboundary = OpenBoundary(grid, South())
eboundary = OpenBoundary(grid, East())
wboundary = OpenBoundary(grid, West())

domain = Subzero.Domain(nboundary, sboundary, eboundary, wboundary)

# Floe instantiation
nfloes = 100
file = jldopen("examples/floe_shapes.jld2", "r")
nfloes = nfloes > size(file["floe_vertices"], 1) ? size(file["floe_vertices"], 1) : nfloes
floe_coords = file["floe_vertices"][1:1]
floe_arr = initialize_floe_field(floe_coords, domain, hmean, Δh)
close(file)

model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
#consts = Constants(E = modulus, Cd_io = 0.0, Cd_ia = 0.0, Cd_ao = 0.0, f = 0.0, μ = 0.0)  # collisions without friction 
simulation = Simulation(model = model, consts = consts, Δt = Δt, nΔt = 4000, COLLISION = true, verbose = true)

# Output setup
initwriter = InitialStateOutputWriter(dir = "output/voronoi", filename = "initial_state.jld2", overwrite = true)
floewriter = FloeOutputWriter(10, dir = "output/voronoi", filename = "f.jld2", overwrite = true)

# Run simulation
run!(simulation, [initwriter, floewriter])