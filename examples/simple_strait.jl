using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero
import LibGEOS as LG

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
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
    two_way_coupling_on = true,
)

# Floe creation
floe_arr = initialize_floe_field(
    FT,
    50,
    [0.7],
    domain,
    hmean,
    Δh,
    rng = Xoshiro(3),
    coupling_settings = coupling_settings,
)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 75,
        npieces = 3,
        nhistory = 1000,
        deform_on = false,
)

# Output setup
dir = "output/simple_strait"
run_time!(simulation) = @time run!(simulation)

initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 3000,
    verbose = false,
    coupling_settings = coupling_settings,
    fracture_settings = fracture_settings,
    writers = writers,
)
    
# Run simulation
run_time!(simulation)

Subzero.create_sim_gif("output/simple_strait/floes.jld2", 
                       "output/simple_strait/initial_state.jld2",
                       "output/simple_strait/simple_strait.gif")