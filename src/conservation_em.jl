using Plots, NCDatasets

jul_path = "output/floe/f.nc"
sim_data = NCDataset(jul_path)
u = sim_data["u"][:, :]
v = sim_data["v"][:, :]
m = sim_data["mass"][:, :]
ξ = sim_data["ξ"][:, :]
moment = sim_data["moment"][:, :]

# Energy conservation
energy = 0.5 * m .* (u.^2 .+ v.^2)
total_energy = sum(energy .+ (0.5 * moment .* ξ.^2), dims = 2)
plot(total_energy, title = "Total Kinetic Energy", legend = false, xlabel = "25 Timesteps", ylabel = "[N]")

# Momentum
momentum = m .* (u .+ v)
total_momentum = sum(momentum, dims = 2)
plot(total_momentum, title = "Total Momentum", legend = false, xlabel = "25 Timesteps", ylabel = "[N * s]")
