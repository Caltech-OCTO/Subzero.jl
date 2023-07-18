using JLD2, Random, Statistics, Subzero, BenchmarkTools, StructArrays
import LibGEOS as LG

Δt = 10
max_overlap = 0.75
Lx = 1e5
Ly = Lx
hmean = 0.25
Δh = 0.0
grid = RegRectilinearGrid(
    (-Lx, Lx),
    (-Ly, Ly),
    1e4,
    1e4,
)
nboundary = PeriodicBoundary(North, grid)
sboundary = PeriodicBoundary(South, grid)
eboundary = CollisionBoundary(East, grid)
wboundary = OpenBoundary(West, grid)
topo_elem = TopographyElement(
    [[[1e4, 0.0], [0.0, 1e4], [1e4, 2e4], [2e4, 1e4], [1e4, 0.0]]],
)
domain = Domain(
    nboundary,
    sboundary,
    eboundary,
    wboundary,
    StructArray([topo_elem]),
)
# Diagonal floe barely overlaping with eastern collision boundary
efloe_small = Floe(
    [[
        [9.5e4, 0.0],
        [9e4, 0.5e4],
        [10e4, 2.5e4],
        [10.05e4, 2e4],
        [9.5e4, 0.0],
    ]],
    hmean,
    Δh,
)
efloe_small.u = 0.5
efloe_small.v = 0.25
consts = Constants()
Subzero.floe_domain_interaction!(efloe_small, domain, consts, Δt, max_overlap)

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 1e4
const hmean = 0.25
const Δh = 0.0
const Δt = 10
# Model instantiation
grid = RegRectilinearGrid(
    (0, Lx),
    (0, Ly),
    Δgrid,
    Δgrid,
)
zero_ocn = Ocean(grid, 0.0, 0.0, 0.0)

zero_atmos = Atmos(grid, 0.0, 0.0, 0.0)


domain = Subzero.Domain(
    CollisionBoundary(North, grid),
    CollisionBoundary(South, grid),
    CompressionBoundary(East, grid, -0.1),
    CompressionBoundary(West, grid, 0.1),
)

# Floe instantiation
coords = [[[0.0, 0.0], [0.0, 1e5], [1e5, 1e5], [1e5, 0.0], [0.0, 0.0]]]

floe_arr = initialize_floe_field(
    FT,
    50,
    [1.0],
    domain,
    0.5,
    0.0,
    min_floe_area = 1e6,
)

model = Model(
    grid,
    zero_ocn,
    zero_atmos,
    domain,
    floe_arr,
)
dir = "output/sim"
writers = OutputWriters(
    InitialStateOutputWriter(
        dir = dir,
        overwrite = true
    ),
    FloeOutputWriter(
        250,
        dir = dir,
        overwrite = true,
    ),
)
simulation = Simulation(
    name = "sim",
    model = model,
    Δt = 10,
    nΔt = 4000,
    writers = writers,
    verbose = true,
    consts = Constants(Cd_io = 0, Cd_ia = 0, Cd_ao = 0),
    fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(model.floes),
        Δt = 75,
    )
)

#@benchmark timestep_sim!(simulation, 10) setup=(sim=deepcopy(simulation))


time_run(simulation) = @time run!(simulation)
time_run(simulation)

Subzero.create_sim_gif("output/sim/floes.jld2", 
                       "output/sim/initial_state.jld2",
                       "output/sim/sim.gif")
# # Run simulation
#time_run(simulation)
#Profile.Allocs.clear()
#@time run!(simulation)
#ProfileView.@profview run!(simulation)
#Profile.Allocs.@profile timestep_sim!(simulation, 1)
# Profile.Allocs.@profile sample_rate=1 Subzero.timestep_coupling!(
#     simulation.model.floes,
#     simulation.model.grid,
#     simulation.model.domain,
#     simulation.model.ocean,
#     simulation.model.atmos,
#     simulation.consts,
#     simulation.coupling_settings,
# )


# Profile.Allocs.@profile sample_rate=1 Subzero.timestep_coupling!(
#     simulation.model,
#     simulation.consts,
#     simulation.coupling_settings,
#     Threads.SpinLock(),
# )

# PProf.Allocs.pprof(from_c = false)
# last(sort(results.allocs, by=x->x.size))
# Subzero.create_sim_gif(
#     joinpath(dir, "floes.jld2"), 
#     joinpath(dir, "initial_state.jld2"),
#     joinpath(dir, "test.gif"),
# )
