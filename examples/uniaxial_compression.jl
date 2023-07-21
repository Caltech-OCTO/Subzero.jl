using JLD2, Random, Statistics, StructArrays, Subzero

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
    (0.0, Lx),
    (0.0, Ly),
    Δgrid,
    Δgrid,
)
ocean = Ocean(grid, 0.0, -0.3, 0.0)
atmos = Atmos(grid, 0.0, 0.0, 0.0)

# Domain creation
nboundary = CompressionBoundary(North, grid, -0.5)
sboundary = CompressionBoundary(South, grid, 0.5)
eboundary = CollisionBoundary(East, grid)
wboundary = CollisionBoundary(West, grid)


domain = Domain(nboundary, sboundary, eboundary, wboundary)

# coupling_settings = CouplingSettings(
#     subfloe_point_generator = SubGridPointsGenerator(grid, 2),
#     two_way_coupling_on = false,
# )

# Floe creation
floe_arr = initialize_floe_field(
    FT,
    50,
    [1.0],
    domain,
    hmean,
    Δh,
    rng = Xoshiro(3),
)

# Model creation
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 75,
        npieces = 3,
        nhistory = 1000,
        deform_on = true,
)

# Output setup
dir = "output/uniaxial_compression"
run_time!(simulation) = @time run!(simulation)

initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
floewriter = FloeOutputWriter(50, dir = dir, overwrite = true)
writers = OutputWriters(initwriter, floewriter)

simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 10000,
    verbose = true,
    fracture_settings = fracture_settings,
    writers = writers,
)
    
# Run simulation
run_time!(simulation)
#ProfileView.@profview run!(simulation)

Subzero.create_sim_gif("output/uniaxial_compression/floes.jld2", 
                       "output/uniaxial_compression/initial_state.jld2",
                       "output/uniaxial_compression/uniaxial_compression.gif")