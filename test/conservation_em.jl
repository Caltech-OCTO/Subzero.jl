#using Plots, NCDatasets, MAT

# mat_path = "/Users/skylargering/src/SubZero_MATLAB/Floes"
# mat_files = readdir(mat_path)
# u_mat = zeros(Float64, 0, 2)
# v_mat = zeros(Float64, 0, 2)
# mass_mat = zeros(Float64, 0, 2)
# ksi_mat = zeros(Float64, 0, 2)
# moment_mat = zeros(Float64, 0, 2)
# force_mat = []
# for file in mat_files
#     if splitext(file)[2] == ".mat"
#         mat_data = matread(joinpath(mat_path, file))
#         u_mat = vcat(u_mat, mat_data["Floe"]["Ui"])
#         v_mat = vcat(v_mat, mat_data["Floe"]["Vi"])
#         mass_mat = vcat(mass_mat, mat_data["Floe"]["mass"])
#         ksi_mat = vcat(ksi_mat, mat_data["Floe"]["ksi_ice"])
#         moment_mat = vcat(moment_mat, mat_data["Floe"]["inertia_moment"])
#         forces = mat_data["Floe"]["collision_force"]
#         force_mat = push!(force_mat, forces)
#     end
# end

# energy_mat = 0.5 * mass_mat .* (u_mat.^2 + v_mat.^2)
# total_energy_mat = sum(energy_mat, dims = 2)
# Δe_mat = (total_energy_mat[end] - total_energy_mat[1])/total_energy_mat[1]

# momentum_mat = mass_mat .* (u_mat .+ v_mat)
# total_momentum_mat = sum(momentum_mat, dims = 2)
# Δm_mat = (total_momentum_mat[end] - total_momentum_mat[1])/total_momentum_mat[1]

function calc_orbital_angular_vel(centroid, u, v)
    x, y = first.(centroid), last.(centroid)
    rad = sqrt.(x.^2 .+ y.^2)
    vel = sqrt.(u.^2 .+ v.^2)
    orbital = zeros(length(x))
    for i in eachindex(x)
        Θ = acos(dot([x[i], y[i]], [u[i], v[i]])/(vel[i]* rad[i]))
        orbital[i] = (vel[i]/rad[i] * sin(Θ))
    end
    return orbital
end

# Energy conservation
function calc_total_energy(u, v, mass, ξ, Ω, moment)
    # Linear kinetic energy
    linear = sum(0.5 * mass .* (u.^2 .+ v.^2))
    # Rotational kinetic energy
    rotational = sum(0.5 * moment .* (ξ.^2 .+ Ω.^2))
    return linear, rotational
end

function calc_total_momentum(u, v, mass, ξ, Ω, moment)
    linear = sum(mass .* (u .+ v))
    angular = sum(moment .* (ξ .+ Ω))
    return linear, angular
end

function check_energy_momentum_conservation(filename, dir)
    file = jldopen(filename, "r")
    tsteps = keys(file["centroid"])
    ntsteps = length(tsteps)
    linear_energy = zeros(ntsteps)
    rotational_energy = zeros(ntsteps)
    linear_momentum = zeros(ntsteps)
    angular_momentum = zeros(ntsteps)
    for i in eachindex(tsteps)
        t = tsteps[i]
        # Needed values
        centroid = file["centroid"][t]
        mass = file["mass"][t]
        moment = file["moment"][t]
        u = file["u"][t]
        v = file["v"][t]
        ξ = file["ξ"][t]
        Ω = calc_orbital_angular_vel(centroid, u, v)
        # calculations
        linear_energy[i], rotational_energy[i] = calc_total_energy(u, v, mass, ξ, Ω, moment)
        linear_momentum[i], angular_momentum[i] = calc_total_momentum(u, v, mass, ξ, Ω, moment)
    end
    close(file)
    # Energy
    total_energy = linear_energy .+ rotational_energy
    plot(
        [total_energy linear_energy rotational_energy],
        title = "Total Kinetic Energy",
        xlabel = "10 Timesteps",
        ylabel = "[N]",
        label=["total" "linear" "rotational"]
    )

    savefig(joinpath(dir, "energy_conservation.png"))
    println("Mean energy: ", mean(total_energy))
    println("Median energy: ", median(total_energy))
    println("Minimum energy: ", minimum(total_energy))
    println("Maximum energy: ", maximum(total_energy))

    # Momentum
    total_momentum = linear_momentum .+ angular_momentum
    plot(
        [total_momentum linear_momentum angular_momentum],
        title = "Total Momentum",
        xlabel = "10 Timesteps",
        ylabel = "[N * s]",
        label=["total" "linear" "angular"]
    )

    savefig(joinpath(dir, "momentum_conservation.png"))
    println("Mean momentum: ", mean(total_momentum))
    println("Median momentum: ", median(total_momentum))
    println("Minimum momentum: ", minimum(total_momentum))
    println("Maximum momentum: ", maximum(total_momentum))
end


using JLD2, Random, SplitApplyCombine, Statistics, StructArrays, Subzero, BenchmarkTools
import LibGEOS as LG

const FT = Float64
const Lx = 1e5
const Ly = Lx
const Δgrid = 1e4
const hmean = 0.25
const Δh = 0.0
const Δt = 10
grid = RegRectilinearGrid(
    FT,
    (0, Lx),
    (0, Ly),
    Δgrid,
    Δgrid,
)

ocean = Ocean(grid, 0.0, 0.0, 0.0)
atmos = Atmos(grid, 0.0, 0.0, 0.0)
nboundary = CollisionBoundary(grid, North())
sboundary = CollisionBoundary(grid, South())
eboundary = CollisionBoundary(grid, East())
wboundary = CollisionBoundary(grid, West())
domain = Domain(nboundary, sboundary, eboundary, wboundary)

rng = Xoshiro(1)
# floe_arr = initialize_floe_field(
#     25,
#     [0.7],
#     domain,
#     hmean,
#     Δh,
#     rng = rng,
#     nhistory = 1000,
#     u = 0.5,
# )
floe_arr = initialize_floe_field(
    [[[[2e4, 2e4], [2e4, 5e4], [5e4, 5e4], [5e4, 2e4], [2e4, 2e4]]],
     [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]],
    domain,
    hmean,
    Δh,
    rng = rng,
    nhistory = 1000,
)
floe_arr.u .= 0.5
#floe_arr = load("output/sim/initial_state.jld2", "sim").model.floes

model = Model(grid, ocean, atmos, domain, floe_arr)

dir = "output/sim"
initwriter = InitialStateOutputWriter(
    dir = dir,
    overwrite = true,
)
floewriter = FloeOutputWriter(
    10,
    dir = dir,
    filename = "floes",
    overwrite = true,
)
writers = OutputWriters(
    initialwriters = StructArray([initwriter]),
    floewriters = StructArray([floewriter]),
)

modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
consts = Constants(E = modulus, μ = 0.0)
coupling_settings = CouplingSettings(coupling_on = false)
fracture_settings = FractureSettings(
        fractures_on = false,
        criteria = HiblerYieldCurve(floe_arr),
        Δt = 75,
        npieces = 3,
        nhistory = 1000,
        deform_on = false,
)
simulation = Simulation(
    model = model,
    consts = consts,
    Δt = Δt,
    nΔt = 1000,
    verbose = true,
    fracture_settings = fracture_settings,
    coupling_settings = coupling_settings,
    writers = writers
)
run!(simulation)
Subzero.create_sim_gif("output/sim/floes.jld2", 
                       "output/sim/initial_state.jld2",
                       "output/sim/test.gif")

check_energy_momentum_conservation("output/sim/floes.jld2", dir)
