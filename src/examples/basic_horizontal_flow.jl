using Subzero, StructArrays
import LibGEOS as LG

# User Inputs
const type = Float64::DataType

const x = 4e5
const y = 4e5
const Δgrid = 10000
const h_mean = 0.5
const Δh = 0.25
const To = 0.0
const Ta = -20.0
const Δt = 10
const newfloe_Δt = 500
ocntemp = 3

# Model instantiation
grid = Grid(x, y, Δgrid)
ocean = Ocean(ones(grid.size), zeros(grid.size), fill(3.0, grid.size))
wind = Wind(ocean)

# Domain creation
domain = Subzero.RectangleDomain(grid, Subzero.PeriodicBC(),
                                     Subzero.PeriodicBC(),
                                     Subzero.CollisionBC(), Subzero.OpenBC())

# Topography instantiation
poly1 = LG.Polygon([[[0.0, 0.0], [5e3, 0.0], [5e3, 1e4],[0.0, 1e4], [0.0, 0.0]]])
topo = Topography(poly1, h_mean)
topo_arr = StructArray(topo for i in 1:1)
# Floe instantiation
poly2 = LG.Polygon([[[0.0, 4e4], [0.0, 1e4], [4e4, 1e4], 
                    [4e4, 4e4], [0.0, 4e4]]])
floe = Floe(poly2, h_mean, Δh)
floe_arr = StructArray(floe for i in 1:1)


model = Model(grid, ocean, wind, domain, topo_arr, floe_arr, To, Ta, Δt, newfloe_Δt)

# Simulation set-up

# show progress

# run!(simulation)

#psi_ocean(xocn, yocn)=0.5e4/1*(sin.(4*(pi/x)*xocn).* sin.(4*(pi/y)*yocn))