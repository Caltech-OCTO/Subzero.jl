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

# Simulation setup
simulation = Simulation(model = model, nΔt = 230)

# Output setup - I think that these should be enums
grid_names = ["u", "v", "du", "dv", "si_frac", "overlap", "mass", "area", "height"]
grid_units = ["m/s", "m/s", "m/s^2", "m/s^2", "unitless", "m^2", "kg", "m^2", "m"]

floe_names = ["centroid", "height", "area", "mass", "moment", "rmax", "coords", "Δα", "u", "v", "ξ", "fxOA", "fyOA", "torqueOA", "p_dxdt", "p_dydt", "p_dudt", "p_dvdt", "p_dξdt", "p_dαdt", "overarea", "alive"]

floe_units = ["location", "m", "m^2", "kg", "kg m^2", "m", "location", "rad", "m/s", "m/s", "rad/s", "N", "N", "N m", "m/s", "m/s", "m/s^2", "m/s^2", "rad/s^2", "rad/s", "m^2", "unitless"]

gridwriter = GridOutputWriter(150, grid_names, grid_units)
floewriter = FloeOutputWriter(100, floe_names, floe_units)
run!(simulation, output_grid, output_data, [gridwriter, floewriter])
