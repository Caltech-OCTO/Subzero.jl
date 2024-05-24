using JLD2, Random, Statistics, Subzero

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(
    (0.0, Lx),
    (0.0, Ly),
    Δgrid,
    Δgrid,
)
ocean = Ocean(grid, 0.0, -0.3, 0.0)
atmos = Atmos(grid, 0.0, 0.0, 0.0)

# Domain creation
nboundary = PeriodicBoundary(North, grid)
sboundary = PeriodicBoundary(South, grid)
eboundary = CollisionBoundary(East, grid)
wboundary = CollisionBoundary(West, grid)

island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]

topo_arr = initialize_topography_field(
    FT,
    [island, topo1, topo2]
)

domain = Domain(nboundary, sboundary, eboundary, wboundary, topo_arr)

coupling_settings = CouplingSettings(
    two_way_coupling_on = true,
)
floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
    stress_calculator = DamageStressCalculator(),
)

# Floe creation
floe_arr = initialize_floe_field(
    FT,
    75,
    [0.7],
    domain,
    hmean,
    Δh,
    Δt,
    rng = Xoshiro(3),
    floe_settings = floe_settings,
)

fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 75,
        npieces = 3,
        deform_on = false,
)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
ridgeraft_settings = RidgeRaftSettings(
    ridge_raft_on = true,
    Δt = 150
)

# Output setup
dir = "output/simple_strait_damage"
run_time!(simulation) = @time run!(simulation)

initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 3000,
    verbose = true,
    floe_settings = floe_settings,
    coupling_settings = coupling_settings,
    fracture_settings = fracture_settings,
    ridgeraft_settings = ridgeraft_settings,
    writers = writers,
)
    
# Run simulation
run_time!(simulation)

plot_sim(
    dir*"/floes.jld2",
    dir*"/initial_state.jld2",
    Δt,
    dir*"/simple_strait.mp4",
)