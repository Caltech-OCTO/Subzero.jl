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
ocean = Ocean(grid, 1.0, 0.0, 3.0)
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
poly2 = LG.Polygon([[[3.2e5, 4e4], [3.2e5, 1e4], [3.7e5, 1e4], 
                    [3.7e5, 4e4], [3.2e5, 4e4]]])
floe1 = Floe(poly2, h_mean, Δh)
floe2 = Floe(Subzero.translate(poly2, [-2e5, 0.0]), h_mean, Δh)
floe_arr = StructArray([floe1, floe2])
consts = Constants()

model = Model(grid, ocean, wind, domain, topo_arr, floe_arr, consts)

output_grid = Grid(model.domain, (10, 8))
output_data = OutputGridData((10, 8))

#plt = Subzero.setup_plot(model)
#cgrid = Grid(model.domain, (coarse_nx, coarse_ny))
#cgrid_data = CoarseGridData(coarse_nx, coarse_ny)

# Simulation set-up
simulation = Simulation(model = model, nΔt = 400)
run!(simulation, output_grid, output_data)