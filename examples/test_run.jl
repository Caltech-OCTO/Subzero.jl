using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, BenchmarkTools, MAT, ProfileView, Profile, PProf
import LibGEOS as LG

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
    FT,
    (-Lx, Lx),
    (-Ly, Ly),
    Δgrid,
    Δgrid,
)
zonal_ocn = Ocean(FT, grid, 0.5, 0.0, 0.0)

zero_atmos = Atmos(grid, 0.0, 0.0, 0.0)

open_domain_no_topo = Subzero.Domain(
    CollisionBoundary(grid, North()),
    CollisionBoundary(grid, South()),
    CollisionBoundary(grid, East()),
    CollisionBoundary(grid, West()),
)

# Floe instantiation

# file = jldopen("test/inputs/floe_shapes.jld2", "r")
# funky_floe_coords = file["floe_vertices"][1:25]
# funky_floe_arr = initialize_floe_field(
#     funky_floe_coords,
#     open_domain_no_topo,
#     hmean,
#     Δh,
# )
# close(file)
funky_floe_arr = initialize_floe_field(
    100,
    [0.5],
    open_domain_no_topo,
    hmean,
    Δh,
    rng = Xoshiro(5),
)
#funky_floe_arr.u .= (-1)^rand(0:1) * (0.1 * rand(length(funky_floe_arr)))
#funky_floe_arr.v .= (-1)^rand(0:1) * (0.1 * rand(length(funky_floe_arr)))

model = Model(
    grid,
    zonal_ocn,
    zero_atmos,
    open_domain_no_topo,
    deepcopy(funky_floe_arr),
)
dir = "output/sim"
writers = OutputWriters(
    # initialwriters = StructArray([InitialStateOutputWriter(
    #     dir = dir,
    #     overwrite = true
    # )]),
    # floewriters = StructArray([FloeOutputWriter(
    #     150,
    #     dir = dir,
    #     overwrite = true,
    # )]),
)
simulation = Simulation(
    name = "sim",
    model = model,
    Δt = 10,
    nΔt = 1500,
    writers = writers,
    verbose = true,
    coupling_settings = CouplingSettings(two_way_coupling_on = false),
    fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(model.floes),
        Δt = 75,
        deform_on = true,
    ),
    rng = Xoshiro(5),
)

#@benchmark timestep_sim!(simulation, 10) setup=(sim=deepcopy(simulation))

# @benchmark Subzero.timestep_collisions!(
#     sim.model.floes,
#     sim.model.max_floe_id,
#     sim.model.domain,
#     zeros(Int, sim.model.max_floe_id),
#     zeros(Int, sim.model.max_floe_id),
#     sim.consts,
#     sim.Δt,
#     sim.collision_settings,
#     Threads.SpinLock(),
# ) setup=(sim=deepcopy(simulation))

# @benchmark Subzero.timestep_coupling!(
#     sim.model,
#     sim.Δt,
#     sim.consts,
#     sim.coupling_settings,
# ) setup=(sim=deepcopy(simulation))

# Subzero.timestep_coupling!(
#     simulation.model,
#     simulation.Δt,
#     simulation.consts,
#     simulation.coupling_settings,
# )

# ProfileView.@profview Subzero.timestep_coupling!(
#     simulation.model.floes,
#     simulation.model.grid,
#     simulation.model.domain,
#     simulation.model.ocean,
#     simulation.model.atmos,
#     simulation.consts,
#     simulation.coupling_settings,
# )

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
