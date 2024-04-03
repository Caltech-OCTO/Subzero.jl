using JLD2, Subzero, Random, Statistics

dir = "output/restart_sim/first_run/"
new_dir = "output/restart_sim/second_run/"

const FT = Float64
const Δt = 10
const nΔt = 5000
const nΔt2 = 5000
const nfloes = 20
const L = 1e5
const Δgrid = 1e4
const hmean = 2
const concentration = 0.7
const uomax = 2

# Run first nΔt steps
grid = RegRectilinearGrid(
    (0.0, L),
    (0.0, L),
    Δgrid,
    Δgrid,
)
ngrid = Int(L/Δgrid) + 1
ygrid = range(0,L,ngrid)

uoprofile = @. uomax * (1 - abs(1 - 2 * ygrid/L))
uvels_ocean = repeat(
    uoprofile,
    outer = (1, ngrid),
)
ocean = Ocean(
    uvels_ocean',
    zeros(grid.Nx + 1, grid.Ny + 1),
    zeros(grid.Nx + 1, grid.Ny + 1),
)

atmos = Atmos(FT, grid, 0.0, 0.0, 0.0)

nboundary = PeriodicBoundary(North, grid)
sboundary = PeriodicBoundary(South, grid)
eboundary = PeriodicBoundary(East, grid)
wboundary = PeriodicBoundary(West, grid)
domain = Domain(nboundary, sboundary, eboundary, wboundary)

floe_settings = FloeSettings(
    subfloe_point_generator = SubGridPointsGenerator(grid, 2),
)
floe_arr = initialize_floe_field(
    FT,
    nfloes,
    [concentration],
    domain,
    hmean,
    0;
    rng = Xoshiro(1),
    floe_settings = floe_settings
)

model = Model(grid, ocean, atmos, domain, floe_arr)

modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus, f = 0, turnθ = 0)

initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
checkpointer = CheckpointOutputWriter(
    250,
    dir = dir,
    filename = "checkpoint.jld2",
    overwrite = true,
    jld2_kw = Dict{Symbol, Any}(),
)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter, checkpointer)

simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = nΔt,
    verbose = true,
    writers = writers,
    rng = Xoshiro(1),
)

run!(simulation)
plot_sim(
    dir * "floes.jld2",
    dir * "initial_state.jld2",
    Δt,
    dir * "shear_flow.mp4",
)

# Run from checkpoint
# is = jldopen(dir * "/initial_state.jld2", "r")
# cp = jldopen(dir * "/checkpoint.jld2",  "r")

# tstep_keys = collect(keys(cp["ocean"]))
# tsteps = unique([if length(filter(isdigit,k)) > 0 parse(Int64, filter(isdigit, k)) end for k in tstep_keys])
# t_max = maximum(tsteps)

# floes = cp["floes"][string(t_max)]
# empty!.(floes.ghosts)

# model2 = Model(
#     is["sim"].model.grid, 
#     cp["ocean"][string(t_max)], 
#     cp["atmos"][string(t_max)], 
#     is["sim"].model.domain, 
#     floes[floes.ghost_id .== 0]
# )

initwriter2 = InitialStateOutputWriter(initwriter; dir = new_dir)
floewriter2 = FloeOutputWriter(floewriter; dir = new_dir)
writers2 = OutputWriters(initwriter2, floewriter2)


Subzero.restart!(dir * "/initial_state.jld2", dir * "/checkpoint.jld2", nΔt2, writers2)
# simulation2 = Simulation(
#     model = model2,
#     consts = is["sim"].consts,
#     Δt = is["sim"].Δt,
#     nΔt = nΔt2,
#     verbose = is["sim"].verbose,
#     writers = writers2,
#     coupling_settings = is["sim"].coupling_settings,
#     simp_settings = is["sim"].simp_settings,
# )

# Subzero.write_data!(simulation2, 0)
# run!(simulation2)
plot_sim(
    new_dir * "floes.jld2",
    new_dir * "initial_state.jld2",
    Δt,
    new_dir * "shear_flow.mp4",
)

# Read new initial_state file and compare size of floes
# is_new = jldopen(new_dir * "/initial_state.jld2", "r")
# size_is = size(is_new["sim"].model.floes)[1]
# size_cp = size(cp["floes"][string(t_max)])[1]

# println("Size of floe vector in checkpoint :: $size_cp")
# println("Size of floe vector in new initial_state AFTER RUNNING SIMULATION :: $size_is")