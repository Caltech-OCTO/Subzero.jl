using Subzero, StructArrays, Statistics, JLD2, SplitApplyCombine
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
grid = RegRectilinearGrid(0, Lx, 0, Ly, Δgrid, Δgrid)
ocean = Ocean(grid, 0.0, -0.2, 0.0)
atmos = Atmos(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), zeros(grid.dims .+ 1))

# Domain creation - boundaries and topography
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]]
#topo = TopographyElement([[[-9.5e4, 4.5e4], [-9.5e4, 6.5e4], [-6.5e4, 6.5e4],
#                           [-6.5e4, 4.5e4], [-9.5e4, 4.5e4]]])
#topo_arr = StructVector([topo for i in 1:1])

topo1 = [[[0, 0.0], [0, 1e5], [2e4, 1e5], [3e4, 5e4], [2e4, 0], [0.0, 0.0]]]
topo2 = [[[8e4, 0], [7e4, 5e4], [8e4, 1e5], [1e5, 1e5], [1e5, 0], [8e4, 0]]]

topo_arr = StructVector([TopographyElement(t) for t in [island, topo1, topo2]])

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe instantiation
#floe1_coords = [[[7.5e4, 7.5e4], [7.5e4, 9.5e4], [9.5e4, 9.5e4], 
#                    [9.5e4, 7.5e4], [7.5e4, 7.5e4]]]
#floe2_coords = [[[6.5e4, 4.5e4], [6.5e4, 6.5e4], [8.5e4, 6.5e4], 
#                 [8.5e4, 4.5e4], [6.5e4, 4.5e4]]]
#floe_arr = StructArray([Floe(c, h_mean, Δh) for c in [floe1_coords, floe2_coords]])
#floe_arr.u[1] = 1
#floe_arr.v[1] = 1
#floe_arr.u[2] = 1
floe_arr = initialize_floe_field(50, [1], domain, 0.5, 0.0)
#model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
#floe1_coords = [[[9.75e4, 7e4], [9.75e4, 5e4], [10.05e4, 5e4], 
#                    [10.05e4, 7e4], [9.75e4, 7e4]]]
#floe2_coords = [[[6.5e4, 4.5e4], [6.5e4, 6.5e4], [8.5e4, 6.5e4], 
#                 [8.5e4, 4.5e4], [6.5e4, 4.5e4]]]
#floe_info = load("output/voronoi_escape_floe.jld2")
#floe_coords = floe_info["coords"]
#floe_mcx = floe_info["mc_x"]
#floe_mcy = floe_info["mc_y"]
# floe_arr = Vector{Floe}()
# for i in eachindex(floe_coords)
#    floe = Floe(floe_coords[i], h_mean, Δh)
#     floe.mc_x = floe_mcx[i]
#     floe.mc_y = floe_mcy[i]
#     push!(floe_arr, floe)
# end
#floe_arr = StructArray([Floe(c, h_mean, Δh) for c in floe_coords])
#floe_arr.mc_x .= floe_mcx
#floe_arr.mc_y .= floe_mcy
#floe_arr.u[1] = 1
#floe_arr.v[1] = 0
#floe_arr.u[2] = 1

#coords = load("output/voronoi/voronoi_escape_floe.jld2")["floes/2246"].coords[1]
#floe = Floe(coords, h_mean, Δh)
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup

#modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants()
#consts = Constants(E = modulus, Cd_io = 0.0, Cd_ia = 0.0, Cd_ao = 0.0, f = 0.0, μ = 0.0)  # collisions without friction 
simulation = Simulation(model = model, consts = consts, Δt = Δt, nΔt = 100, COLLISION = true, verbose = true)

# Output setup
initwriter = InitialStateOutputWriter(dir = "output/voronoi", filename = "initial_state.jld2", overwrite = true)
gridwriter = GridOutputWriter(50, grid, (10, 10), dir = "output/voronoi", filename = "g.nc", overwrite = true)
floewriter = FloeOutputWriter(10, dir = "output/voronoi", filename = "f.jld2", overwrite = true)
checkpointwriter = CheckpointOutputWriter(2246, dir = "output/voronoi", overwrite = true)

# Run simulation
run!(simulation, [initwriter, floewriter])