using JLD2, Random, Statistics, Subzero

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.125
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(
    (0.0, Lx),
    (0.0, Ly),
    Δgrid,
    Δgrid,
)
ocean = Ocean(grid, 0.1, 0.0, 0.0)
atmos = Atmos(grid, 0.0, 0.0, -1.0)

# Domain creation
nboundary = CompressionBoundary(North, grid, -0.1)
sboundary = CompressionBoundary(South, grid, 0.1)
eboundary = CollisionBoundary(East, grid)
wboundary = CollisionBoundary(West, grid)

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe creation
floe_arr = initialize_floe_field(
    FT,
    75,
    [1.0],
    domain,
    hmean,
    Δh,
    rng = Xoshiro(1),
)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
run_time!(simulation) = @time run!(simulation)
dir = "output/packed_compression"

# Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(40, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

# Simulation settings 
ridgeraft_settings = RidgeRaftSettings(
    ridge_raft_on = true,
    Δt = 150,
    domain_gain_probability = 0.5
)
coupling_settings = CouplingSettings(two_way_coupling_on = true)


simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 20000,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    coupling_settings = coupling_settings,
    ridgeraft_settings = ridgeraft_settings,
)
run_time!(simulation)

Subzero.plot_sim(
    "output/packed_compression/floes.jld2",
    "output/packed_compression/initial_state.jld2",
    20,
    "output/packed_compression/packed_compression.mp4",
)