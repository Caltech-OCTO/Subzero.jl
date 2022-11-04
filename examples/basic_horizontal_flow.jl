using Subzero, StructArrays
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
ocean = Ocean(grid, 0.5, 0.5, 0.0)
wind = Wind(zeros(grid.dims), zeros(grid.dims), fill(-20.0, grid.dims))

# Domain creation
domain = Subzero.RectangleDomain(grid, northBC = OpenBC(),
                                       southBC = OpenBC(),
                                       eastBC = OpenBC(),
                                       westBC = OpenBC())

                                    
# Topography instantiation - TODO - figure out how to not have topography...
poly1 = LG.Polygon([[[0.0, -0.75e5], [5e3, -0.75e5], [5e3, -1e5],[0.0, -1e5], [0.0, -0.75e5]]])
topo = Topography(poly1, h_mean)
topo_arr = StructVector([topo for i in 1:1])
# Floe instantiation
floe_poly = LG.Polygon([[[-0.5e4, 1e4], [-0.5e4, -1e4], [0.5e4, -1e4], 
                    [0.5e4, 1e4], [-0.5e4, 1e4]]])
floe1 = Floe(floe_poly, h_mean, Δh)
floe_arr = StructArray([floe1 for i in 1:1])
consts = Constants()

model = Model(grid, ocean, wind, domain, topo_arr, floe_arr, consts)

# Simulation setup
simulation = Simulation(model = model, nΔt = 7200)

# Output setup
gridwriter = GridOutputWriter([GridOutput(i) for i in 1:9], 10, "g.nc", grid, (10, 10))
floewriter = FloeOutputWriter([FloeOutput(i) for i in 1:24], 10, "nwf.nc", grid)

# Run simulation
run!(simulation, [])

#TODO: We should add dissolved mass
