using Subzero, StructArrays, Statistics, JLD2, SplitApplyCombine
import LibGEOS as LG

# User Inputs
const type = Float64::DataType
#const Lx = 1e5
#const Ly = 1e5
const Δgrid = 10000
const h_mean = 0.25
const Δh = 0.0
const Δt = 20
const newfloe_Δt = 500
const coarse_nx = 10
const coarse_ny = 10

Lx = 1e5
Ly = 1e5
grid = RegRectilinearGrid(-Lx, Lx, -Lx, Lx, 1e4, 1e4)
double_periodic_domain = Domain(PeriodicBoundary(grid, North()), PeriodicBoundary(grid, South()),
                                PeriodicBoundary(grid, East()), PeriodicBoundary(grid, West()))
coords1 = splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))
coords2 = splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], -[7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))
floe_arr = StructArray(Floe([c], 0.5, 0.0) for c in [coords1, coords2])
for i in eachindex(floe_arr)
    floe_arr.id[i] = i
end
trans_arr = StructArray([Floe(Subzero.translate([coords1], [0.0, -2Ly]), 0.5, 0.0),
                            Floe(Subzero.translate([coords2], [2Lx, 0.0]), 0.5, 0.0)])
for i in eachindex(trans_arr)
    trans_arr.id[i] = i
end
Subzero.timestep_collisions!(trans_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
add_ghosts!(floe_arr, double_periodic_domain)
Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
# Model instantiation
grid = RegRectilinearGrid(-Lx, Lx, 0, Ly, Δgrid, Δgrid)
ocean = Ocean(grid, 1.0, 0.0, 0.0)
atmos = Atmos(zeros(grid.dims .+ 1), zeros(grid.dims .+ 1), fill(0.0, grid.dims .+ 1))

# Domain creation - boundaries and topography
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())

topo = TopographyElement([[[-9.5e4, 4.5e4], [-9.5e4, 6.5e4], [-6.5e4, 6.5e4],
                           [-6.5e4, 4.5e4], [-9.5e4, 4.5e4]]])
topo_arr = StructVector([topo for i in 1:1])

domain = Domain(nboundary, sboundary, eboundary, wboundary)

# Floe instantiation
floe1_coords = [[[9.75e4, 7e4], [9.75e4, 5e4], [10.05e4, 5e4], 
                    [10.05e4, 7e4], [9.75e4, 7e4]]]
#floe2_coords = [[[6.5e4, 4.5e4], [6.5e4, 6.5e4], [8.5e4, 6.5e4], 
#                 [8.5e4, 4.5e4], [6.5e4, 4.5e4]]]
floe_arr = StructArray([Floe(c, h_mean, Δh) for c in [floe1_coords]])
floe_arr.u[1] = 1
floe_arr.v[1] = 0
#floe_arr.u[2] = 1
model = Model(grid, ocean, atmos, domain, floe_arr)

# Simulation setup
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
#consts = Constants(E = modulus, Cd_io = 0.0, Cd_ia = 0.0, Cd_ao = 0.0, f = 0.0, μ = 0.0)  # collisions without friction 
simulation = Simulation(model = model, consts = consts, Δt = Δt, nΔt = 3000, COLLISION = true)

# Output setup
gridwriter = GridOutputWriter([GridOutput(i) for i in 1:9], 10, "g.nc", grid, (10, 10))
floewriter = FloeOutputWriter([FloeOutput(i) for i in [3:4; 6; 9:11; 22:26]], 30, "f.nc", grid)

# Run simulation
run!(simulation, [floewriter])