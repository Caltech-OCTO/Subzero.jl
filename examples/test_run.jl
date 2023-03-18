using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, BenchmarkTools, MAT
import LibGEOS as LG

# User Inputs
const FT = Float64
const Lx = 5e4
const Ly = 7e4
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const Δt = 10

# Model instantiation
grid = RegRectilinearGrid(
    FT,
    (-1e4, Lx),
    (0, Ly),
    Δgrid,
    Δgrid,
)
ocean = Ocean(grid, 0.0, 0.0, 0.0)
atmos = Atmos(grid, 0.0, 0.0, 0.0)


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
floe1 = [[[2e4, 2e4], [2e4, 5e4], [5e4, 5e4], [5e4, 2e4], [2e4, 2e4]]]
floe_on_wall = initialize_floe_field(
        [Subzero.translate(floe1, [-2.9e4, 0.0])],
        domain, # Just affects shape, type doesnt' matter
        0.25,
        0.0,
        rng = rng,
        nhistory = 100,
)
floe_on_wall.u[1] = -0.1
floe_on_wall.v[1] = 0.01
# file = jldopen("examples/floe_shapes.jld2", "r")
#nfloes = nfloes > size(file["floe_vertices"], 1) ? size(file["floe_vertices"], 1) : nfloes
# floe_coords = file["floe_vertices"][3:5]
# floe_coords[3] = Subzero.translate(floe_coords[3], [0.0, 2e4])
# floe_arr = initialize_floe_field(floe_coords, domain, hmean, Δh, rng = rng)
# close(file)
# floe_arr.u[1] = 0.1
# floe_arr.v[2] = -0.2
# floe_arr.v[3] = 0.2
#floe_arr.u .= rand(rng, nfloes) * 0.05
#floe_arr.v .= rand(rng, nfloes) * 0.05
#floe_arr.ξ .= rand(rng, nfloes) * 1e-12

# floe_arr = initialize_floe_field(
#     10,
#     [0.5],
#     domain,
#     hmean,
#     Δh,
#     rng = rng,
#     nhistory = 1000,
# )
# nfloes = length(floe_arr)
# floe_arr.u .= rand(rng, nfloes) * 0.1
# floe_arr.v .= rand(rng, nfloes) * 0.1
#floe_arr.ξ .= rand(rng, nfloes) * 1e-6
# floe_arr = initialize_floe_field(
#     [[[[2e4, 2e4], [2e4, 5e4], [5e4, 5e4], [5e4, 2e4], [2e4, 2e4]]],
#     [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]],
#     [[[5.5e4, 2e4], [5.25e4, 4e4], [5.75e4, 4e4], [5.5e4, 2e4]]]],
#     domain,
#     hmean,
#     Δh,
#     rng = rng,
#     nhistory = 1000,
# )
# floe_arr.u[1] = 0.11
# floe_arr.u[2] = -0.1
# floe_arr.ξ[3] = 1e-5
# floe_arr.v[1] = 0.001
# floe_arr.v[2] = 0.001
# floe_arr.v[3] = 0.001

#[[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]
#floe_arr = load("output/sim/thread1_initial_state.jld2", "sim").model.floes

model = Model(grid, ocean, atmos, domain, floe_on_wall)

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
    nΔt = 10000,
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

em_lists = Subzero.check_energy_momentum_conservation_julia(
    joinpath(dir, string(tstring, "_floes.jld2")),
    dir,
)
println(
    [vals[1] != 0 ?
        (vals[end] - vals[1])/vals[1] * 100 :
        @warn "Starting value is 0"
    for vals in em_lists]
)

Subzero.create_sim_gif(
    joinpath(dir, string(tstring, "_floes.jld2")), 
    joinpath(dir, string(tstring, "_initial_state.jld2")),
    joinpath(dir, string(tstring, "_test.gif")),
)
