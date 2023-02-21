using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero
import LibGEOS as LG

# User Inputs
const type = Float64::DataType
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
grid = RegRectilinearGrid(0, Lx, 0, Ly, Δgrid, Δgrid)
uvels = repeat(hcat(range(0.1 ,0.6, step = 0.1)', range(0.5, 0.1, step = -0.1)'), outer = 11)
ocean = Ocean(uvels, zeros(grid.dims .+ 1), zeros(grid.dims .+ 1))
atmos = Atmos(grid, 0.0, 0.0, -1.0)

# Domain creation
nboundary = PeriodicBoundary(grid, North())
sboundary = PeriodicBoundary(grid, South())
eboundary = PeriodicBoundary(grid, East())
wboundary = PeriodicBoundary(grid, West())

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe creation
floe_arr = initialize_floe_field(50, [0.8], domain, 0.25, 0.0, rng = Xoshiro(1))

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
simulation = Simulation(model = model, consts = consts, Δt = Δt, nΔt = 7000, COLLISION = true, verbose = true)

# Output setup
dir = "output/converge_diverge"
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
checkpointwriter = CheckpointOutputWriter(1000, dir = dir, overwrite = true)

# Run simulation
run!(simulation, [initwriter, floewriter, checkpointwriter])

Subzero.create_sim_gif("output/converge_diverge/floes.jld2", 
                       "output/converge_diverge/initial_state.jld2",
                       "output/converge_diverge/sheer_flow.gif")