using Subzero, Random, Statistics

const FT = Float64
const hmean = 0.25
const Δh = 0.0
const Δy = 500.0
const Δx = 500.0
Nx = 256
Ny = 256
const Lx = Δx*Nx
const Ly = Δy*Ny 

save_time_interval = 60
save_time_interval_checkpoint = 3600*24*10
stop_time = 60*200
Δt = 60

# Model instantiation
floe_grid = Subzero.RegRectilinearGrid(
    FT,
    (0, Lx),
    (0, Ly),
    Δx,
    Δx,
)

nboundary = PeriodicBoundary(North, floe_grid)
sboundary = PeriodicBoundary(South, floe_grid)
eboundary = PeriodicBoundary(East, floe_grid)
wboundary = PeriodicBoundary(West, floe_grid)
domain = Subzero.Domain(nboundary, sboundary, eboundary, wboundary)

dir = "output/test_writers"
save_it_interval = save_time_interval / Δt
initwriter = InitialStateOutputWriter(dir = dir, filename = "initial_state.jld2", overwrite = true)
floewriter = FloeOutputWriter(save_it_interval, dir = dir, filename = "floes.jld2", overwrite = true)
checkpointwriter = Subzero.CheckpointOutputWriter(save_time_interval_checkpoint, dir = dir, filename = "checkpoint.jld2", overwrite = true)
writers = Subzero.OutputWriters(initwriter, floewriter, checkpointwriter)

ocean = Subzero.Ocean(FT,floe_grid, 0.0, 0.0, 0.0)
atmos = Subzero.Atmos(floe_grid, 0.0, 0.0, -1.0)
floe_arr = initialize_floe_field(50, [0.40], domain, hmean, Δh, rng = Xoshiro(1))

floe_model = Subzero.Model(floe_grid, ocean, atmos, domain, floe_arr)

modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus)
floe_simulation = Subzero.Simulation(
    model = floe_model,
    consts = consts,
    Δt = Δt,
    nΔt = Int(floor(stop_time/Δt)),
    verbose = true,
    writers = writers,
    coupling_settings = CouplingSettings(two_way_coupling_on = true),
)
tstep = 0
while tstep <= floe_simulation.nΔt
    Subzero.timestep_sim!(floe_simulation, tstep)
    global tstep += 1
end

f = Figure()
Axis(f[1, 1])
h = heatmap!(vec(xgrid), vec(ygrid), vec(simulation.model.ocean.τx'), colormap = :rainbow)
for i in eachindex(simulation.model.floes)
    x = cos(simulation.model.floes.α[i]).*simulation.model.floes.mc_x[i] .- sin(simulation.model.floes.α[i]).*simulation.model.floes.mc_y[i]
    y = sin(simulation.model.floes.α[i]).*simulation.model.floes.mc_x[i] .+ cos(simulation.model.floes.α[i]).*simulation.model.floes.mc_y[i]
    x .+= simulation.model.floes.centroid[i][1]
    y .+= simulation.model.floes.centroid[i][2]
    scatter!(x, y)
end
Colorbar(f[1, 2], h)
f