using Subzero, StructArrays
import LibGEOS as LG

# User Inputs
const type = Float64::DataType

const Lx = 4e5
const Ly = 3e5
const Δgrid = 10000
const h_mean = 0.5
const Δh = 0.25
const Δt = 10
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

# Model instantiation
grid = Grid(Lx, Ly, Δgrid, Δgrid)
ocean = Ocean(grid, 2.0, 0.0, 3.0)
wind = Wind(ones(grid.dims), zeros(grid.dims), fill(4.0, grid.dims))

# Domain creation
domain = Subzero.RectangleDomain(grid, northBC = OpenBC(),
                                       southBC = OpenBC(),
                                       eastBC = OpenBC(),
                                       westBC = OpenBC())

                                    
# Topography instantiation
poly1 = LG.Polygon([[[0.0, 0.0], [5e3, 0.0], [5e3, 1e4],[0.0, 1e4], [0.0, 0.0]]])
topo = Topography(poly1, h_mean)
topo_arr = StructArray(topo for i in 1:1)

# Floe instantiation
poly2 = LG.Polygon([[[3.45e5, 4e4], [3.45e5, 1e4], [3.95e5, 1e4], 
                    [3.95e5, 4e4], [3.45e5, 4e4]]])
poly3 = LG.Polygon([[[2e5, 4e4], [2.5e5, 4e4], [2.5e5, 3e4], [2e5, 4e4]]])
floe1 = Floe(poly2, h_mean, Δh)
floe2 = Floe(poly3, h_mean, Δh)
floe_arr = StructArray([floe1, floe2])
consts = Constants()

model = Model(grid, ocean, wind, domain, topo_arr, floe_arr, consts)

# Simulation setup
simulation = Simulation(model = model, nΔt = 500)

# Output setup
gridwriter = GridOutputWriter(vcat([GridOutput(i) for i in 1:4], GridOutput(9)), 100, "g.nc", grid, (9, 10))
floewriter = FloeOutputWriter([GridOutput(i) for i in 1:10], 100, "f.nc", grid)

# Run simulation
run!(simulation, [gridwriter, floewriter])
