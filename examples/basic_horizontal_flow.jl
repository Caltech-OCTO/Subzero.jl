using Subzero, StructArrays
import LibGEOS as LG




floe_poly4 = LG.Polygon([[[2., 2.], [8., 2.], [8., 8.], [2., 2.]]])
floe4 = Subzero.Floe(floe_poly4, 0.25, 0.0)
area_ratio4, _, _, idx1 = Subzero.floe_area_ratio(floe4, collect(-5.:5.:10.), collect(0.:5.:10.))


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
domain = Subzero.RectangleDomain(grid, northBC = PeriodicBC(),
                                       southBC = PeriodicBC(),
                                       eastBC = CollisionBC(),
                                       westBC = OpenBC())

                                    
# Topography instantiation
poly1 = LG.Polygon([[[0.0, 0.0], [5e3, 0.0], [5e3, 1e4],[0.0, 1e4], [0.0, 0.0]]])
topo = Topography(poly1, h_mean)
topo_arr = StructArray(topo for i in 1:1)

# Floe instantiation
poly2 = LG.Polygon([[[0.0, 4e4], [0.0, 1e4], [4e4, 1e4], 
                    [4e4, 4e4], [0.0, 4e4]]])
floe = Floe(poly2, h_mean, Δh)
floe_arr = StructArray(floe for i in 1:1)


model = Model(grid, ocean, wind, domain, topo_arr, floe_arr, Δt, newfloe_Δt)

#plt = Subzero.setup_plot(model)
#cgrid = Grid(model.domain, (coarse_nx, coarse_ny))
#cgrid_data = CoarseGridData(coarse_nx, coarse_ny)

# Simulation set-up
simulation = Simulation(model = model)
run!(simulation)

Subzero.floe_area_ratio(floe, grid.xg, grid.yg)