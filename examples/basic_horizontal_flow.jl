using Subzero, StructArrays, Statistics
import LibGEOS as LG

# User Inputs
const type = Float64::DataType

const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const h_mean = 0.25
const Δh = 0.0
const Δt = 10
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

# Model instantiation
grid = Grid(-Lx, Lx, -Ly, Ly, Δgrid, Δgrid)
ocean = Ocean(grid, 1.0, 0.0, 0.0)
wind = Wind(zeros(grid.dims), zeros(grid.dims), fill(-20.0, grid.dims))

# Domain creation
domain = Subzero.RectangleDomain(grid, northBC = CollisionBC(),
                                       southBC = CollisionBC(),
                                       eastBC = CollisionBC(),
                                       westBC = CollisionBC())
                                    
# Topography instantiation - TODO - figure out how to not have topography...
poly1 = LG.Polygon([[[0.0, -0.75e5], [5e3, -0.75e5], [5e3, -1e5],[0.0, -1e5], [0.0, -0.75e5]]])
topo = Topography(poly1, h_mean)
topo_arr = StructVector([topo for i in 1:1])
# Floe instantiation
floe1_poly = LG.Polygon([[[1.25e4, 9e4], [1.25e4, -9e4], [3.75e4, -9e4], 
                    [3.75e4, 9e4], [1.25e4, 9e4]]])
floe1 = Floe(floe1_poly, h_mean, Δh, u = 0.0)
floe2_poly = LG.Polygon([[[0.7e4, 8e4], [0.7e4, 6e4], [1.2e4, 6e4], 
                    [1.2e4, 8e4], [0.7e4, 8e4]]])
floe2 = Floe(floe2_poly, h_mean, Δh, u = 0.0)
floe_arr = StructArray([floe2])

modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
#consts = Constants(E = modulus, Cd_io = 0.0, Cd_ia = 0.0, Cd_ao = 0.0, f = 0.0, μ = 0.0)  # collisions without friction 
consts = Constants(E = modulus)
model = Model(grid, ocean, wind, domain, topo_arr, floe_arr, consts)

# Simulation setup
simulation = Simulation(model = model, Δt = Δt, nΔt = 3000, COLLISION = true)

# Output setup
gridwriter = GridOutputWriter([GridOutput(i) for i in 1:9], 10, "g.nc", grid, (10, 10))
floewriter = FloeOutputWriter([FloeOutput(i) for i in [3:4; 6; 9:11; 22:26]], 100, "f.nc", grid)

# Run simulation
run!(simulation, [floewriter])