using JLD2, Subzero, Random, Statistics

const FT = Float64
const Δt = 10
const nΔt = 5000
const n_part_sim = 3
const nfloes = 20
const L = 1e5
const Δgrid = 1e4
const hmean = 2
const concentration = 0.7
const uomax = 2

dirs = [joinpath("output/restart_sim", "run_" * string(i)) for i in 1:n_part_sim]

# Build initial model / simulation
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

floe_settings = FloeSettings(subfloe_point_generator = SubGridPointsGenerator(grid, 2))
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

initwriter = InitialStateOutputWriter(dir = dirs[1], overwrite = true)
checkpointer = CheckpointOutputWriter(
    250,
    dir = dirs[1],
    filename = "checkpoint.jld2",
    overwrite = true,
    jld2_kw = Dict{Symbol, Any}(),
)
floewriter = FloeOutputWriter(50, dir = dirs[1], overwrite = true)
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

# Run the first part of the simulation
run!(simulation)

# Run remaining parts of the simulation
for i in 2:n_part_sim
    # Build new output writers
    global initwriter = InitialStateOutputWriter(initwriter; dir = dirs[i])
    global floewriter = FloeOutputWriter(floewriter; dir = dirs[i])
    writers = if i < n_part_sim
        global checkpointer = CheckpointOutputWriter(checkpointer; dir = dirs[i])
        OutputWriters(initwriter, floewriter, checkpointer)
    else
        OutputWriters(initwriter, floewriter)
    end
    # Run next part of the simulation
    Subzero.restart!(
        dirs[i-1] * "/initial_state.jld2",
        dirs[i-1] * "/checkpoint.jld2",
        nΔt, writers; start_tstep = nΔt * (i - 1))
end

# Plot all simulation parts
for i in 1:n_part_sim
    plot_sim(
        joinpath(dirs[i], "floes.jld2"),
        joinpath(dirs[i], "initial_state.jld2"),
        Δt,
        joinpath(dirs[i], "shear_flow.mp4"),
    )
    println("Simulation part $i plotted.")
end