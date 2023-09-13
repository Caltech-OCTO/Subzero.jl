using JLD2, Random, Statistics, Subzero, BenchmarkTools, StructArrays, SplitApplyCombine
import LibGEOS as LG

grid = RegRectilinearGrid(
           (0, 2e5),
           (0, 2e5),
           1e4,
           1e4,
       )
domain = Subzero.Domain(
           CollisionBoundary(North, grid),
           CollisionBoundary(South, grid),
           CollisionBoundary(East, grid),
           CollisionBoundary(West, grid),
)
consts = Constants()
coords = [[[[0.1e4, 0.1e4], [0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [0.1e4, 0.1e4]]], [[[1.9e4, 1.9e4], [1.9e4, 4e4], [4e4, 4e4], [4e4, 1.9e4], [1.9e4, 1.9e4]]]]
floes = initialize_floe_field(
                  Float64,
                  coords,
                  domain,
                  1.0,
                  0.0,
)
Subzero.timestep_collisions!(
    floes,
    2,
    domain,
    consts,
    10,
    CollisionSettings(floe_floe_max_overlap = 0.99),
    Threads.SpinLock(),
)
Subzero.timestep_ridging_rafting!(floes, 2, domain,
    Subzero.RidgeRaftSettings(ridge_probability = 1.0, raft_probability = 0.0),
    Subzero.CouplingSettings(),
    SimplificationSettings(),
    consts,
    10,
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
