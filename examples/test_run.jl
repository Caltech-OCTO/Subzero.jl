using JLD2, Random, Statistics, Subzero, BenchmarkTools, StructArrays, SplitApplyCombine, Test
import LibGEOS as LG

A = [[[39990.01810326673, 91800.1335126005], [40255.74085266778, 92422.63261641243], [42260.0005801607, 96996.25896636098], [50276.1908284482, 91234.71709074081], [50376.42934215918, 86715.43624628006], [49362.76249022356, 85145.19992670821], [48664.61560665339, 84082.84702554355], [45900.518277471376, 84621.4307064748], [44689.21871417973, 85381.91995212912], [40280.55275261199, 89878.91285449876], [40262.40862608565, 89869.5115308193], [40244.76588360351, 89887.41477207805], [40224.73555929154, 89876.97024145136], [40122.36902960179, 89980.13143162905], [40120.268228183246, 89979.02673121012], [40034.469411763836, 90065.39285187308], [39990.01810326673, 91800.1335126005]]]
B = [ [[27976.81255056931, 82050.49102583545], [29414.04069428343, 84244.99448129853], [33839.104850523334, 86538.67992233086], [39704.6052143005, 89579.31454598375], [39704.58292620962, 89580.47585011902], [40280.55275261199, 89878.91285449876], [44898.78399147662, 85168.1560669079], [35114.94447362642, 80776.13880559478], [33226.32322119895, 81045.64614300567], [27976.81255056931, 82050.49102583545]]]
Subzero.find_shared_edges_midpoint(A, B)

frac_floe = deepcopy(frac_deform_floe)  # Without interactions, won't deform
no_frac_floe = Floe(  # This floe is colliding with frac_deform_floe
    [[
        [1467.795, -25319.563],
        [1664.270, -25640.216],
        [-1105.179, -33458.936],
        [-17529.019, -50035.583],
        [-21193.828, -50088.777],
        [-21370.170, -32618.322],
        [-21247.656, -31077.536],
        [-12818.593, -27031.048],
        [1467.795, -25319.563],
    ]],
    0.25,
    0.0,
)
no_frac_small = Floe(  # This floe is too small to fracture or deform
    [[
        [1e3, 1e3],
        [1e3, 1.5e3],
        [1.5e3, 1.5e3],
        [1.5e3, 1e3],
        [1e3, 1e3],
    ]],
    0.25,
    0.0,
)
frac_deform_floe.stress = frac_stress
frac_deform_floe.interactions = collect([
    3,
    -279441968.984,
    -54223517.438,
    -21091.0918258529,
    -40358.0042297616,
    -148920620521.112,
    6795329.38154967,
]')
frac_deform_floe.num_inters = 1
frac_deform_floe.p_dudt = 0.11
frac_floe.stress = frac_stress
no_frac_small.stress = frac_stress

floes = StructArray([
    frac_deform_floe, frac_floe, no_frac_floe, no_frac_small
])
floes.id .= collect(1:4)
new_floes = Subzero.split_floe(
    floes[1],
    Xoshiro(3),
    FractureSettings(
        fractures_on = true,
        npieces = 2,
        criteria = HiblerYieldCurve(floes),
        Δt = 75,
        deform_on = true,
    ),
    CouplingSettings(),
    Constants(),
    10,
) 

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
# bounds_overlap_area = LibGEOS.area(LibGEOS.intersection(
#     LibGEOS.Polygon(floes_base.coords[1]),
#     boundary_poly,
# ))
# topo_overlap_area = LibGEOS.area(LibGEOS.intersection(
#     LibGEOS.Polygon(floes_base.coords[2]),
#     topo_poly,
# ))
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
rafting_settings = Subzero.RidgeRaftSettings(
    ridge_probability = 0.0,  # force ridging
    raft_probability = 1.0,
)
# first floe overlaps with both other floes, and it will break when ridged with floe 2

coords = [
        [[[-0.1e4, -0.1e4], [-0.1e4, 2e4], [2e4, 2e4], [2e4, -0.1e4], [-0.1e4, -0.1e4]]],
        [[[3e4, 3e4], [3e4, 5e4], [5e4, 5e4], [5e4, 3e4], [3e4, 3e4]]]
    ]
floes_base = setup_floes_with_inters(coords, collision_domain, consts,
    collision_settings, lock,
)
floes = deepcopy(floes_base)
update_height(floes, 1, 0.1, consts)  # floe 1 will ridge onto floe 2
update_height(floes, 2, 0.1, consts)
total_mass = floes.mass[1] + floes.mass[2]
h1, h2 = floes.height
area1, area2 = floes.area
cent1, cent2 = floes.centroid
pieces_list = StructArray{Floe{Float64}}(undef, 0)

max_id = Subzero.timestep_ridging_rafting!(
    floes,
    pieces_list,
    2,
    collision_domain,
    maximum(floes.id),
    rafting_settings,
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
