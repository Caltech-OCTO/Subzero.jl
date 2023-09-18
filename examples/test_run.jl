using JLD2, Random, Statistics, Subzero, BenchmarkTools, StructArrays, SplitApplyCombine, Test
import LibGEOS as LG

function setup_floes_with_inters(coords, domain, consts,
    collision_settings, lock,  Δx = nothing, Δy = nothing,
)
    floes = initialize_floe_field(
        Float64,
        coords,
        domain,
        1.0,
        0.0,
    )
    if !isnothing(Δx)
        for i in eachindex(Δx)
            Subzero.translate!(floes.coords[i], Δx[i], Δy[i])
            floes.centroid[i][1] += Δx[i]
            floes.centroid[i][2] += Δy[i]
        end
    end
    add_ghosts!(floes, domain)
    Subzero.timestep_collisions!(  # Add interactions
        floes,
        length(floes),
        domain,
        consts,
        10,
        collision_settings, 
        lock,
    )
    return floes
end
function update_height(floes, i, new_height, consts)
    floes.height[i] = new_height
    floes.mass[i] = floes.area[i] * floes.height[i] * consts.ρi
    floes.moment[i] = Subzero.calc_moment_inertia(
        floes.coords[i],
        floes.centroid[i],
        floes.height[i],
        ρi = consts.ρi,
    )
end

grid = RegRectilinearGrid(
    (0, 1e5),
    (0, 1e5),
    1e4,
    1e4,
)
topo_coords = [[[5e4, 5e4], [5e4, 7e4], [7e4, 7e4], [7e4, 5e4], [5e4, 5e4]]]
collision_domain = Subzero.Domain(
    CollisionBoundary(North, grid),
    CollisionBoundary(South, grid),
    CollisionBoundary(East, grid),
    CollisionBoundary(West, grid),
    initialize_topography_field([topo_coords])
)
periodic_domain = Subzero.Domain(
    PeriodicBoundary(North, grid),
    PeriodicBoundary(South, grid),
    PeriodicBoundary(East, grid),
    PeriodicBoundary(West, grid),
)

consts = Constants()
coupling_settings = CouplingSettings()
simp_settings = SimplificationSettings()
collision_settings = CollisionSettings(floe_floe_max_overlap = 0.99) # don't fuse
lock = Threads.SpinLock()

ridge_settings = Subzero.RidgeRaftSettings(
    ridge_probability = 1.0,  # no ridging
    raft_probability = 0.0,  # no rafting
)
# first floe overlaps with both other floes, and it will break when ridged with floe 2

coords = [
    [[[8e4, 0.75e4], [10e4, 2.75e4], [10.5e4, 2.75e4], [8.5e4, 0.75e4], [8e4, 0.75e4]]],
    [[[9e4, 0.1e4], [9e4, 2.25e4], [9.5e4, 2.25e4], [9.5e4, 0.1e4], [9e4, 0.1e4]]],
    [[[0.1e4, 0.1e4], [0.1e4, 2.25e4], [2.25e4, 2.25e4], [2.25e4, 0.1e4], [0.1e4, 0.1e4]]],
]
base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
    collision_settings, lock
)
floes = deepcopy(base_floes)
update_height(floes, 1, 0.1, consts)
update_height(floes, 4, 0.1, consts)  # ghost of floe 1
total_mass = sum(floes.mass[1:3])
h1, h2, h3, h4 = floes.height
area1, area2, area3, area4 = floes.area
cent1, cent2, cent3, cent4 = floes.centroid
pieces_list = StructArray{Floe{Float64}}(undef, 0)
max_id = Subzero.timestep_ridging_rafting!(
    floes,
    pieces_list,
    3,
    collision_domain,
    maximum(floes.id),
    ridge_settings,
    coupling_settings,
    simp_settings,
    consts,
    10,
)

# Ridging with domain
ridge_settings = Subzero.RidgeRaftSettings(
    ridge_probability = 1.0,  # force ridging
    raft_probability = 0.0,
)
update_height(floes, 1, 1.0, consts)  # floe1 will ridge onto floe 2
total_mass = sum(floes.mass)
Subzero.timestep_ridging_rafting!(
    floes,
    1,
    domain,
    ridge_settings,
    coupling_settings,
    simp_settings,
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
