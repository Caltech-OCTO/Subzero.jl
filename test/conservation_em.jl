using Plots, NCDatasets, MAT

mat_path = "/Users/skylargering/src/SubZero_MATLAB/Floes"
mat_files = readdir(mat_path)
u_mat = zeros(Float64, 0, 2)
v_mat = zeros(Float64, 0, 2)
mass_mat = zeros(Float64, 0, 2)
ksi_mat = zeros(Float64, 0, 2)
moment_mat = zeros(Float64, 0, 2)
force_mat = []
for file in mat_files
    if splitext(file)[2] == ".mat"
        mat_data = matread(joinpath(mat_path, file))
        u_mat = vcat(u_mat, mat_data["Floe"]["Ui"])
        v_mat = vcat(v_mat, mat_data["Floe"]["Vi"])
        mass_mat = vcat(mass_mat, mat_data["Floe"]["mass"])
        ksi_mat = vcat(ksi_mat, mat_data["Floe"]["ksi_ice"])
        moment_mat = vcat(moment_mat, mat_data["Floe"]["inertia_moment"])
        forces = mat_data["Floe"]["collision_force"]
        force_mat = push!(force_mat, forces)
    end
end

energy_mat = 0.5 * mass_mat .* (u_mat.^2 + v_mat.^2)
total_energy_mat = sum(energy_mat, dims = 2)
Δe_mat = (total_energy_mat[end] - total_energy_mat[1])/total_energy_mat[1]

momentum_mat = mass_mat .* (u_mat .+ v_mat)
total_momentum_mat = sum(momentum_mat, dims = 2)
Δm_mat = (total_momentum_mat[end] - total_momentum_mat[1])/total_momentum_mat[1]

jul_path = "output/floe/f.nc"
sim_data = NCDataset(jul_path)
u = sim_data["u"][:, :]
v = sim_data["v"][:, :]
m = sim_data["mass"][:, :]
ξ = sim_data["ξ"][:, :]
moment = sim_data["moment"][:, :]
xforce = sim_data["xcollision_force"][:, :]
yforce = sim_data["ycollision_force"][:, :]


# Energy conservation
energy = 0.5 * m .* (u.^2 .+ v.^2)
total_energy = sum(energy, dims = 2)
plot(total_energy, title = "Total Kinetic Energy", legend = false, xlabel = "25 Timesteps", ylabel = "[N]")
Δe = (total_energy[end] - total_energy[1])/total_energy[1]

# Momentum
momentum = m .* (u .+ v)
total_momentum = sum(momentum, dims = 2)
plot(total_momentum, title = "Total Momentum", legend = false, xlabel = "25 Timesteps", ylabel = "[N * s]")
Δm = (total_momentum[end] - total_momentum[1])/total_momentum[1]

# Make gif
Subzero.create_sim_gif(jul_path, collect(-1e5:1e4:1e5), collect(-1e5:1e4:1e5))