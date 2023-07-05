using JLD2, Random, Statistics, Subzero, BenchmarkTools, StructArrays
import LibGEOS as LG

frac_stress = [-29955.396 -3428.008; -3428.008	-1942.0464]
frac_deform_floe = Floe(
    [[
        [-50548.186, -49995.968],
        [-50550.745, -37790.078],
        [-20856.010, -32518.566],
        [-20929.577, -49989.757],
        [-50548.186, -49995.968],
    ]],
    0.25,
    0.0,
    u = 0.1,
    v = -0.2,
    ξ = 0.05,
)
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
frac_deform_floe.p_dudt = 0.11
frac_floe.stress = frac_stress
no_frac_small.stress = frac_stress

floes = StructArray([
    frac_deform_floe, frac_floe, no_frac_floe, no_frac_small
])
floes.id .= collect(1:4)
frac_settings = FractureSettings(
    fractures_on = true,
    criteria = HiblerYieldCurve(floes),
    Δt = 75,
    deform_on = true,
)
new_floes = Subzero.split_floe(
    floes[1],
    Xoshiro(3),
    frac_settings,
    CouplingSettings(),
    Constants(),
    10,
) 

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
