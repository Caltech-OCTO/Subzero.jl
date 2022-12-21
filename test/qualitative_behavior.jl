"""
Simulations that show simple, fundamental behavior of sea ice.
Run and check these simulations to confirm behavior has not changed. 
"""

using Subzero, StructArrays

const T = Float64::DataType
const Lx = 1e5
const Ly = 1e5
const Δgrid = 10000
const h_mean = 0.25
const Δh = 0.0
const nΔt = 4000
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10
const floe_fn = "f.nc"

grid = RegRectilinearGrid(-Lx, Lx, -Ly, Ly, Δgrid, Δgrid)

"""
Simulation:
Expected Behavior:
"""
zero_ocn = Ocean(grid, 0.0, 0.0, 0.0)
meridional_ocn = Ocean(grid, 0.0, 1.0, 0.0)

zero_wind = Wind(grid, 0.0, 0.0, 0.0)
zonal_wind = Wind(grid, 3.0, 0.0, 0.0)

open_domain_no_topo = Subzero.Domain(OpenBoundary(grid, North()),
                                     OpenBoundary(grid, South()),
                                     OpenBoundary(grid, East()), 
                                     OpenBoundary(grid, West()))

topography = TopographyElement([[[-3e4, 0.0], [-3e4, 2e4], [-3.5e4, 2e4], 
                                 [-3.5e4, 0.0], [-3e4, 0.0]]])
collision_domain_topo = Subzero.Domain(CollisionBoundary(grid, North()),
                                       CollisionBoundary(grid, South()),
                                       CollisionBoundary(grid, East()), 
                                       CollisionBoundary(grid, West()),
                                       StructArray([topography]))

stationary_rect_floe = StructArray([Floe([[[0.0, 0.0], [0.0, 2e4], [0.5e4, 2e4],
                                          [0.5e4, 0.0], [0.0, 0.0]]], h_mean, Δh)])

zonal_2rect_floe = StructArray([Floe([[[-5e4, 0.0], [-5e4, 2e4], [-5.5e4, 2e4],
                                          [-5.5e4, 0.0], [-5e4, 0.0]]], h_mean, Δh, u = 3.0),
                                Floe([[[0.0, 0.0], [0.0, 2e4], [0.5e4, 2e4],
                                          [0.5e4, 0.0], [0.0, 0.0]]], h_mean, Δh, u = -3.0),
                                Floe([[[3.e4, 0.0], [3e4, 2e4], [3.5e4, 2e4],
                                          [3.5e4, 0.0], [3e4, 0.0]]], h_mean, Δh, u = 0.0)])
standard_consts = Constants()
no_ocndrag_consts = Constants(Cd_io = 0.0, Cd_ao = 0.0, μ = 0.0)
no_drag_consts = Constants(Cd_io = 0.0, Cd_ao = 0.0, Cd_ia = 0.0, μ = 0.0, f = 0.0)

floewriter = FloeOutputWriter([FloeOutput(i) for i in [3:4; 6; 9:11; 22:26]], 30, floe_fn, grid)

sim_arr = Vector{Simulation}(undef, 0)
"""
Simulation 1: One floe pushed by meridional (south-to-north) 1m/s ocean flow. Floe is initally stationary.
Expected Behavior: Floe velocity quickly reaches ocean velocity and flows northward. 
"""
model1 = Model(grid, meridional_ocn, zero_wind, open_domain_no_topo, deepcopy(stationary_rect_floe))
simulation1 = Simulation(name = "sim1", model = model1, consts = standard_consts, Δt = 10, nΔt = nΔt, COLLISION = false)
push!(sim_arr, simulation1)
"""
Simulation 2: One floe pushed by zonal (west-to-east) 1m/s wind flow. Ocean-ice drag coefficent set to 0.
            Floe is initally stationary.
Expected Behavior: Floe should drift to the right due to Coriolis force in the northern hemisphere. 
"""
model2 = Model(grid, zero_ocn, zonal_wind, open_domain_no_topo, deepcopy(stationary_rect_floe))
simulation2 = Simulation(name = "sim2", model = model2, consts = no_ocndrag_consts, Δt = 50, nΔt = nΔt, COLLISION = false)
push!(sim_arr, simulation2)

"""
Simulation 3: One floe with initial velocity flows into a collision boundary and then bounces off a topography element.
Another floe woth opposite initial velocity hits the other side of the collision boundary. No drag so the floes can move
free of the ocean.
Expected Behavior: Floes should bounce off of the topography elements and then off of the walls. 
"""
model3 = Model(grid, zero_ocn, zero_wind, collision_domain_topo, deepcopy(zonal_2rect_floe))
simulation3 = Simulation(name = "sim3", model = model3, consts = no_drag_consts, Δt = 10, nΔt = nΔt, COLLISION = true)
push!(sim_arr, simulation3)

# Run the simulations
for sim in sim_arr
    run!(sim, [floewriter])
    Subzero.create_sim_gif(joinpath("output", sim.name, floe_fn),
                           joinpath("output", sim.name, "domain.jld2"),
                           joinpath("output", sim.name, string(sim.name, ".gif")))
end