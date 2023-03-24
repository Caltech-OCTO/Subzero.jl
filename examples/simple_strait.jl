using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero
import LibGEOS as LG

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(
    FT,
    (0, Lx),
    (0, Ly),
    Δgrid,
    Δgrid,
)
ocean = Ocean(grid, 0.0, -0.3, 0.0)
atmos = Atmos(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), zeros(grid.dims .+ 1))

# Domain creation
nboundary = PeriodicBoundary(grid, North())
sboundary = PeriodicBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]

topo_arr = StructVector([TopographyElement(t) for t in [island, topo1, topo2]])
domain = Domain(nboundary, sboundary, eboundary, wboundary, topo_arr)

# Floe creation
floe_arr = initialize_floe_field(50, [0.7], domain, hmean, Δh, rng = Xoshiro(1))

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
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
checkpointwriter = CheckpointOutputWriter(1000, dir = dir, overwrite = true)
writers = OutputWriters(
    initialwriters = StructArray([initwriter]),
    floewriters = StructArray([floewriter]),
    checkpointwriters = StructArray([checkpointwriter]),
)

simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 3000,
    verbose = true,
    fracture_settings = fracture_settings,
    writers = writers,
)

# Run simulation
run!(simulation)

Subzero.create_sim_gif("output/simple_strait/floes.jld2", 
                       "output/simple_strait/initial_state.jld2",
                       "output/simple_strait/simple_strait.gif")