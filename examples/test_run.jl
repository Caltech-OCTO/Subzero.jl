using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, BenchmarkTools, MAT
import LibGEOS as LG

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const Δt = 10

# Model instantiation
grid = RegRectilinearGrid(
    FT,
    (0, Lx),
    (0, Ly),
    Δgrid,
    Δgrid,
)
ocean = Ocean(grid, 0.25, 0.0, 0.0)
atmos = Atmos(
    zeros(grid.dims .+ 1),
    zeros(grid.dims .+ 1),
    zeros(grid.dims .+ 1),
)

# Domain creation - boundaries and topography
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
topo = TopographyElement([[[-9.5e4, 4.5e4], [-9.5e4, 6.5e4], [-6.5e4, 6.5e4],
                          [-6.5e4, 4.5e4], [-9.5e4, 4.5e4]]])
topo_arr = StructVector([topo for i in 1:1])

topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]

topo_arr = StructVector([TopographyElement(t) for t in [island, topo1, topo2]])

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe instantiation

rng = Xoshiro(1)
floe_arr = initialize_floe_field(
    10,
    [0.5],
    domain,
    hmean,
    Δh,
    rng = rng,
    nhistory = 1000,
)
nfloes = length(floe_arr)
floe_arr.u .= rand(rng, nfloes) * 0.1
floe_arr.v .= rand(rng, nfloes) * 0.1
#floe_arr.ξ .= rand(rng, nfloes) * 1e-6
# floe_arr = initialize_floe_field(
#     [[[[2e4, 3e4], [2e4, 6e4], [5e4, 6e4], [5e4, 3e4], [2e4, 3e4]]],
#      [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]],
#     domain,
#     hmean,
#     Δh,
#     rng = rng,
#     nhistory = 1000,
# )
# floe_arr.u[1] = 0.1
# floe_arr.u[2] = -0.1
#floe_arr = load("output/sim/thread1_initial_state.jld2", "sim").model.floes

model = Model(grid, ocean, atmos, domain, floe_arr)

# Output setup
#tstring = string("thread", Threads.nthreads())
tstring = "offset_blocks"
#tstring = "thread1_rerun"
dir = "output/sim"
initwriter = InitialStateOutputWriter(
    dir = dir,
    filename = string(tstring, "_initial_state.jld2"),
    overwrite = true,
)
# gridwriter = GridOutputWriter(
#     100,
#     model.grid,
#     (10, 10),
#     dir = dir,
#     filename = string(tstring, "_grid.nc"),
#     overwrite = true,
# )
floewriter = FloeOutputWriter(
    10,
    dir = dir,
    filename = string(tstring, "_floes.jld2"),
    overwrite = true,
)
# checkpointwriter = CheckpointOutputWriter(
#     500,
#     dir = dir,
#     filename = string(tstring, "_checkpointer.jld2"),
#     overwrite = true,
# )

writers = OutputWriters(
    initialwriters = StructArray([initwriter]),
    floewriters = StructArray([floewriter]),
)
# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus, μ = 0.0)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 5000,
    verbose = true,
    fracture_settings = FractureSettings(
        fractures_on = false,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 75,
        npieces = 3,
        nhistory = 1000,
        deform_on = false,
    ),
    coupling_settings = CouplingSettings(
        coupling_on = false,
    ),
    writers = writers,
)


time_run(simulation) = @time run!(simulation)

# Run simulation
time_run(simulation)

Subzero.check_energy_momentum_conservation_julia(joinpath(dir, string(tstring, "_floes.jld2")), dir)

Subzero.create_sim_gif(
    joinpath(dir, string(tstring, "_floes.jld2")), 
    joinpath(dir, string(tstring, "_initial_state.jld2")),
    joinpath(dir, string(tstring, "_test.gif")),
)
