using Subzero, StructArrays, Statistics
import LibGEOS as LG

# User Inputs
const type = Float64::DataType

const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const h_mean = 0.25
const Δh = 0.0
const Δt = 20
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

# Model instantiation
grid = Grid(-Lx, Lx, -Ly, Ly, Δgrid, Δgrid)
ocean = Ocean(grid, 0.0, 0.0, 0.0)
wind = Wind(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), fill(-20.0, grid.dims .+ 1))

# Domain creation - boundaries and topography
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

topo = TopographyElement([[[0.5e4, 5e4], [0.5e4, 7e4], [1e4, 7e4], [1e4, 5e4], [0.5e4, 5e4]]])
topo_arr = StructVector([topo for i in 1:1])

domain = Subzero.Domain(nboundary, sboundary, eboundary, wboundary, topo_arr)

# Floe instantiation
floe1_poly = LG.Polygon([[[-1.75e4, 5e4], [-1.75e4, 7e4], [-1.25e4, 7e4], 
                    [-1.25e4, 5e4], [-1.75e4, 5e4]]])
floe1 = Floe(floe1_poly, h_mean, Δh, u = 1.0)
floe_arr = StructArray([floe1])

model = Model(grid, ocean, wind, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
#consts = Constants(E = modulus)
consts = Constants(E = modulus, Cd_io = 0.0, Cd_ia = 0.0, Cd_ao = 0.0, f = 0.0, μ = 0.0)  # collisions without friction 
simulation = Simulation(model = model, consts = consts, Δt = Δt, nΔt = 6000, COLLISION = true)

# Output setup
gridwriter = GridOutputWriter([GridOutput(i) for i in 1:9], 10, "g.nc", grid, (10, 10))
floewriter = FloeOutputWriter([FloeOutput(i) for i in [3:4; 6; 9:11; 22:26]], 30, "f.nc", grid)

# Run simulation
run!(simulation, [floewriter])