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
# Set u velocities
uvels = zeros(FT, (Nx + 1, Ny + 1))
for i in CartesianIndices(uvels)
    r, c = Tuple(i)
    if c ≤ 2
        uvels[i] = 0.2
    elseif c ≥ Nx
        uvels[i] = -0.2
    end
end
# Set v velocities
vvels = zeros(FT, (Nx + 1, Ny + 1))
for i in CartesianIndices(vvels)
    r, c = Tuple(i)
    if r ≤ 2
        vvels[i] = -0.2
    elseif r ≥ Nx
        vvels[i] = 0.2
    end
end

ocean = Ocean(
    uvels,
    vvels,
    zeros(grid.Nx + 1, grid.Ny + 1),
)
atmos = Atmos(grid, 0.0, 0.0, -1.0)

# Domain creation
nboundary = OpenBoundary(North, grid)
sboundary = OpenBoundary(South, grid)
eboundary = OpenBoundary(East, grid)
wboundary = OpenBoundary(West, grid)

domain = Domain(nboundary, sboundary, eboundary, wboundary)

floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
)
# Floe creation
floe_arr = initialize_floe_field(
    FT,
    50,
    [0.8],
    domain,
    hmean,
    Δh;
    rng = Xoshiro(1),
    floe_settings = floe_settings
)

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
    nΔt = 10000,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
    coupling_settings = coupling_settings,
)
run_time!(simulation)

plot_sim(
    "output/shear_flow/floes.jld2",
    "output/shear_flow/initial_state.jld2",
    20,
    "output/shear_flow/shear_flow.mp4",
)