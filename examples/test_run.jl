using JLD2, Random, Statistics, Subzero, BenchmarkTools, StructArrays, SplitApplyCombine
import LibGEOS as LG


Δt = 10
Lx = 1e5
Ly = 1e5
collision_settings = CollisionSettings()
grid = RegRectilinearGrid(
    (-Lx, Lx),
    (-Lx, Lx),
    1e4,
    1e4,
)
double_periodic_domain = Domain(
    PeriodicBoundary(North, grid),
    PeriodicBoundary(South, grid),
    PeriodicBoundary(East, grid),
    PeriodicBoundary(West, grid),
)
# Parent-parent collison (parents are touching)
coords1 = splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))
coords2 = splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], -[7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))

floe_arr = StructArray(Floe([c], 0.5, 0.0) for c in [coords1, coords2])
for i in eachindex(floe_arr)
    floe_arr.id[i] = Float64(i)
end
trans_arr = StructArray([
    Floe(
        Subzero.translate([coords1],
        0.0, -2Ly),
        0.5,
        0.0,
    ),
    Floe(
        Subzero.translate([coords2], 2Lx, 0.0),
        0.5,
        0.0,
    ),
])
for i in eachindex(trans_arr)
    trans_arr.id[i] = i
end
spinlock = Threads.SpinLock()
Subzero.timestep_collisions!(
    trans_arr,
    2,
    double_periodic_domain,
    Subzero.Constants(),
    Δt,
    collision_settings,
    spinlock,
)
add_ghosts!(floe_arr, double_periodic_domain)
Subzero.timestep_collisions!(
    floe_arr,
    2,
    double_periodic_domain,
    Subzero.Constants(),
    Δt,
    collision_settings,
    spinlock,
)


# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
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
zonal_ocn = Ocean(grid, 0.5, 0.0, 0.0)

zero_atmos = Atmos(grid, 0.5, 0.0, 0.0)


domain = Subzero.Domain(
    PeriodicBoundary(North, grid),
    PeriodicBoundary(South, grid),
    PeriodicBoundary(East, grid),
    PeriodicBoundary(West, grid),
)

# Floe instantiation
coords = [[
    [0.0, 0.0],
    [5e4, 1e5],
    [1e5, 1e5],
    [1e5, 0.0],
    [0.0, 0.0],
]]

coupling_settings = CouplingSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 1),
    two_way_coupling_on = true,
)

floe_arr = initialize_floe_field(
    FT,
    30,
    [0.8],
    domain,
    0.5,
    0.0;
    coupling_settings = coupling_settings,
    rng = Xoshiro(1),
)
floe_arr.α .= 1e-5
model = Model(
    grid,
    zonal_ocn,
    zero_atmos,
    domain,
    floe_arr,
)
dir = "output/sim"
writers = OutputWriters(
    # InitialStateOutputWriter(
    #     dir = dir,
    #     overwrite = true
    # ),
    # FloeOutputWriter(
    #     250,
    #     dir = dir,
    #     overwrite = true,
    # ),
)
simulation = Simulation(
    name = "sim",
    model = model,
    Δt = 10,
    nΔt = 200,
    writers = writers,
    verbose = true,
    consts = Constants(Cd_ia = 0, Cd_ao = 0),
    coupling_settings = coupling_settings,
    
)

#@benchmark timestep_sim!(simulation, 10) setup=(sim=deepcopy(simulation))
function run_with_saving(simulation, tstep_between_save)
    nsave = fld(simulation.nΔt, tstep_between_save) + 1
    τx = zeros(FT, grid.Nx + 1, grid.Ny + 1, nsave)
    Subzero.startup_sim(simulation, nothing, 1)
    tstep = 0
    counter = 1
    while tstep <= simulation.nΔt
        # Timestep the simulation forward
        timestep_sim!(simulation, tstep)
        if mod(tstep, tstep_between_save) == 0
            τx[:, :, counter] .= simulation.model.ocean.τx
            counter += 1
        end
        tstep+=1
    end
    Subzero.teardown_sim(simulation)
    return τx
end

τx = run_with_saving(simulation, 50)


# time_run(simulation) = @time run!(simulation)
# time_run(simulation)

# Subzero.create_sim_gif("output/sim/floes.jld2", 
#                        "output/sim/initial_state.jld2",
#                        "output/sim/sim.gif")
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
