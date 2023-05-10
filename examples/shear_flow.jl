using JLD2, Random, Statistics, StructArrays, Subzero
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
ocean = Ocean(FT, uvels, zeros(grid.dims .+ 1), zeros(grid.dims .+ 1))
atmos = Atmos(grid, 0.0, 0.0, -1.0)

# Domain creation
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe creation
floe_arr = initialize_floe_field(50, [0.8], domain, hmean, Δh, rng = Xoshiro(1))

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)

# Run simulation
run_time!(simulation) = @time run!(simulation)
dir = "output/shear_flow"
for i in 1:10
    # Output setup
    local initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
    local floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)

    local writers = OutputWriters(
        initialwriters = StructArray([initwriter]),
        floewriters = StructArray([floewriter]),
    )

    local simulation = Simulation(
        model = model,
        consts = consts,
        Δt = Δt,
        nΔt = 4320,
        verbose = false,
        writers = writers,
        rng = Xoshiro(1),
    )
    run_time!(simulation)
end
 
Subzero.create_sim_gif("output/shear_flow/floes.jld2", 
                       "output/shear_flow/initial_state.jld2",
                       "output/shear_flow/shear_flow.gif")