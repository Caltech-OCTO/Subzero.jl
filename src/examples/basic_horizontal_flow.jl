using Subzero
import LibGEOS as LG

# User Inputs
const type = Float64::DataType

const x = 4e5
const y=4e5
const dgrid = 10000
const h_mean = 0.5
const h_delta = 0.25
const icetemp = -20.0
const ocntemp = 3.0
const dt = 10
const newfloe_dt = 500

# Model instantiation
grid = Grid(x, y, dgrid)
ocean = Ocean(ones(grid.size), zeros(grid.size), fill(ocntemp, grid.size))
wind = Wind(ocean)
poly = LG.Polygon([[[0.0, 4e4], [0.0, 1e4], [4e4, 1e4], 
                    [4e4, 4e4], [0.0, 4e4]]])
floe = Floe(poly, h_mean, h_delta)
#floe_arr = Subzero.StructArrays.StructArray([floe])
# Simulation set-up
mean_heatflux, iceheight = Subzero.calc_icestuff(ocean, -3, 10, 500)

model = Model(grid, ocean, wind, floe, mean_heatflux, iceheight, 920.0, 1.4e-4, 15*pi/180)
# show progress

# run!(simulation)

#psi_ocean(xocn, yocn)=0.5e4/1*(sin.(4*(pi/x)*xocn).* sin.(4*(pi/y)*yocn))