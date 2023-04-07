using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, JLD2
import LibGEOS as LG

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(
    FT,
    (0, Lx),
    (0, Ly),
    Δgrid,
    Δgrid,
)
uvels = repeat([range(0, 0.5, length = 26); range(0.5, 0, length = 25)], outer = (1, 51))
ocean = Ocean(uvels, zeros(grid.dims .+ 1), zeros(grid.dims .+ 1))
atmos = Atmos(grid, 0.0, 0.0, -1.0)

# Domain creation
nboundary = PeriodicBoundary(grid, North())
sboundary = PeriodicBoundary(grid, South())
eboundary = PeriodicBoundary(grid, East())
wboundary = PeriodicBoundary(grid, West())

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe creation
function Base.convert(::Type{Subzero.StressCircularBuffer{Float64}}, cb::DataStructures.CircularBuffer{Matrix{Float64}})
    return Subzero.StressCircularBuffer{Float64}(cb, sum(cb))
end
file = jldopen("output/sheer_flow/golden_initial_state.jld2", "r")
floe_coords = file["sim"].model.floes.coords
floe_mcx = file["sim"].model.floes.mc_x
floe_mcy = file["sim"].model.floes.mc_y
close(file)
floe_arr = initialize_floe_field(floe_coords, domain, hmean, Δh)
floe_arr.mc_x .= floe_mcx
floe_arr.mc_y .= floe_mcy
#floe_arr = initialize_floe_field(50, [0.8], domain, hmean, Δh, rng = Xoshiro(1))

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Output setup
dir = "output/sheer_flow"
initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
checkpointwriter = CheckpointOutputWriter(1000, dir = dir, overwrite = true)

writers = OutputWriters(
    initialwriters = StructArray([initwriter]),
    floewriters = StructArray([floewriter]),
    checkpointwriters = StructArray([checkpointwriter]),
)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 7000,
    verbose = true,
    writers = writers,
)

# Run simulation
run!(simulation)

Subzero.create_sim_gif("output/sheer_flow/floes.jld2", 
                       "output/sheer_flow/initial_state.jld2",
                       "output/sheer_flow/sheer_flow.gif")