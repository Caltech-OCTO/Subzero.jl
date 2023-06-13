using JLD2, Random, Statistics, StructArrays, Subzero
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
uvels = repeat(
    [range(0, 0.5, length = 26); range(0.5, 0, length = 25)],
    outer = (1, 51),
)
ocean = Ocean(
    uvels',
    zeros(grid.Nx + 1, grid.Ny + 1),
    zeros(grid.Nx + 1, grid.Ny + 1),
)
atmos = Atmos(grid, 0.0, 0.0, -1.0)

# Domain creation
nboundary = PeriodicBoundary(North, grid)
sboundary = PeriodicBoundary(South, grid)
eboundary = PeriodicBoundary(East, grid)
wboundary = PeriodicBoundary(West, grid)

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe creation
floe_arr = initialize_floe_field(FT, 50, [0.8], domain, hmean, Δh, rng = Xoshiro(1))

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
run_time!(simulation) =  @time run!(simulation)
dir = "output/shear_flow"

    # Output setup
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 2000,
    verbose = false,
    writers = writers,
    rng = Xoshiro(1),
)
run_time!(simulation)
 
# Subzero.create_sim_gif("output/shear_flow/floes.jld2", 
#                        "output/shear_flow/initial_state.jld2",
#                        "output/shear_flow/shear_flow.gif")