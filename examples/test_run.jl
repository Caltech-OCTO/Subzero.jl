using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, BenchmarkTools, MAT
import LibGEOS as LG

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 1e4
const hmean = 0.25
const Δh = 0.0
const Δt = 10

grid = RegRectilinearGrid(
    Float64,
    (-2.5e4, 1e5),
    (-2.5e4, 1e5),
    1e4,
    1e4,
)
open_domain_no_topo = Subzero.Domain(
    OpenBoundary(FT, North, grid),
    OpenBoundary(FT, South, grid),
    OpenBoundary(FT, East, grid),
    OpenBoundary(FT, West, grid),
)
coords1 = [[  # large floe
    [0.0, 0.0],
    [0.0, 1e4],
    [1e4, 1e4],
    [1e4, 0.0],
    [0.0, 0.0],
]]
coords2 = [[  # small floe
    [8e3, 5e3],
    [8e3, 8e3],
    [1.2e4, 8e3],
    [1.2e4, 5e3],
    [8e3, 5e3],
]]
coords3 = [[  # large floe
    [1.1e4, 0.0],
    [1.1e4, 1e4],
    [2.1e4, 1e4],
    [2.1e4, 0.0],
    [1.1e4, 0.0]
]]
coords4 = [[  # small floe
    [5e3, -2e3],
    [5e3, 3e3],
    [8e3, 3e3],
    [8e3, -2e3],
    [5e3, -2e3],
]]

floe_arr = initialize_floe_field(
    FT,
    [coords1, coords2, coords3, coords4],
    open_domain_no_topo,
    0.5,
    0.0,
    min_floe_area = 1e6,
    rng = Xoshiro(1),
)
floe_arr.status[1].tag = Subzero.fuse
floe_arr.status[1].fuse_idx = [2, 4]
floe_arr.status[2].tag = Subzero.fuse
floe_arr.status[2].fuse_idx = [1, 3]
floe_arr.status[3].tag = Subzero.fuse
floe_arr.status[3].fuse_idx = [2]
floe_arr.status[4].tag = Subzero.fuse
floe_arr.status[4].fuse_idx = [4]

max_floe_id = Subzero.fuse_floes!(
    floe_arr,
    4,
    CouplingSettings(),
    10,
    Constants(),
    Xoshiro(1),
)





# Model instantiation
grid = RegRectilinearGrid(
    FT,
    (-Lx, Lx),
    (-Ly, Ly),
    Δgrid,
    Δgrid,
)
zonal_ocn = Ocean(FT, grid, 0.5, 0.0, 0.0)

zero_atmos = Atmos(FT, grid, 0.0, 0.0, 0.0)


domain = Subzero.Domain(
    CollisionBoundary(FT, North, grid),
    CollisionBoundary(FT, South, grid),
    CollisionBoundary(FT, East, grid),
    CollisionBoundary(FT, West, grid),
)

# Floe instantiation
f = jldopen("examples/floe_shapes.jld2", "r") 
funky_floe_arr = initialize_floe_field(
    FT,
    vec(f["floe_vertices"]),
    domain,
    0.25,
    0.0
)
close(f)

Subzero.simplify_floes!(
    funky_floe_arr,
    domain.topography,
    SimplificationSettings(tol = 100.0),
    CollisionSettings(),
    CouplingSettings(),
    Constants(),
    Xoshiro(5),
)
# funky_floe_arr = initialize_floe_field(
#     FT,
#     100,
#     [0.5],
#     domain,
#     hmean,
#     Δh,
#     rng = Xoshiro(5),
# )
#funky_floe_arr.u .= (-1)^rand(0:1) * (0.1 * rand(length(funky_floe_arr)))
#funky_floe_arr.v .= (-1)^rand(0:1) * (0.1 * rand(length(funky_floe_arr)))

model = Model(
    grid,
    zonal_ocn,
    zero_atmos,
    domain,
    funky_floe_arr,
)
dir = "output/sim"
writers = OutputWriters(
    initialwriters = StructArray([InitialStateOutputWriter(
        dir = dir,
        overwrite = true
    )]),
    floewriters = StructArray([FloeOutputWriter(
        250,
        dir = dir,
        overwrite = true,
    )]),
)
simulation = Simulation(
    name = "sim",
    model = model,
    Δt = 10,
    nΔt = 2000,
    writers = writers,
    verbose = true,
    fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(model.floes),
        Δt = 75,
    )
)

#@benchmark timestep_sim!(simulation, 10) setup=(sim=deepcopy(simulation))

# @benchmark Subzero.timestep_collisions!(
#     sim.model.floes,
#     sim.model.max_floe_id,
#     sim.model.domain,
#     sim.consts,
#     sim.Δt,
#     sim.collision_settings,
#     Threads.SpinLock(),
# ) setup=(sim=deepcopy(simulation))

time_run(simulation) = @time run!(simulation)
time_run(simulation)
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
