using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, BenchmarkTools
import LibGEOS as LG

# User Inputs
const type = Float64::DataType
const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const hmean = 0.25
const Δh = 0.0
const Δt = 20

# Model instantiation
grid = RegRectilinearGrid(0, Lx, 0, Ly, Δgrid, Δgrid)
ocean = Ocean(grid, 0.0, 0.0, 0.0)
atmos = Atmos(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), zeros(grid.dims .+ 1))

# Domain creation - boundaries and topography
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

#island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
#topo = TopographyElement([[[-9.5e4, 4.5e4], [-9.5e4, 6.5e4], [-6.5e4, 6.5e4],
#                           [-6.5e4, 4.5e4], [-9.5e4, 4.5e4]]])
#topo_arr = StructVector([topo for i in 1:1])

#topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
#topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]

#topo_arr = StructVector([TopographyElement(t) for t in [island, topo1, topo2]])

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe instantiation
floe_arr = initialize_floe_field(50, [1], domain, 0.5, 0.0, rng = Xoshiro(1), nhistory = 1000)
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup

modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 50,
    verbose = true,
    fracture_settings = FractureSettings(
        fractures_on = true,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 1,
        npieces = 3,
        nhistory = 1000,
        deform_on = true,
    ),
    coupling_settings = CouplingSettings(
        coupling_on = true,
        calc_ocnτ_on = false,
    )
)

# Output setup
initwriter = InitialStateOutputWriter(dir = "output/sim", filename = "initial_state.jld2", overwrite = true)
gridwriter = GridOutputWriter(50, grid, (10, 10), dir = "output/sim", filename = "grid.nc", overwrite = true)
floewriter = FloeOutputWriter(50, dir = "output/sim", filename = "floes.jld2", overwrite = true)
checkpointwriter = CheckpointOutputWriter(50, dir = "output/sim", overwrite = true)


@benchmark timestep_sim!(simulation, 0, [], Float64)
# Run simulation
#run!(simulation, [floewriter, initwriter])

#Subzero.create_sim_gif("output/sim/floes.jld2", 
#                       "output/sim/initial_state.jld2",
#                       "output/sim/test.gif")