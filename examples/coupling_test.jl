using Subzero
using Statistics
using JLD2

path = "output/coupling_tests"

const FT = Float64;
# Create environment
grid = RegRectilinearGrid(
    FT,
  (0, 1e5), # xbounds
  (0, 1e5), # ybounds
  2e3,      # grid cell width
  2e3,      # grid cell height
 ) 

#  val = LinRange(1.2, -1.2,11);
#  arr = repeat(val, 11);
#  arr2 = transpose(reshape(arr, (11,11))); # converging ocean field

uvels = repeat(
    [range(0, 0.48, length = 25); 0.5; range(0.48, 0, length = 25)],
    outer = (1, 51),
)

ocean = Ocean(
    FT,
    uvels', #converging u velocity
    #fill(0.0, 11, 11),
    fill(0.0, 51, 51), # 0.0m/s v-velocity
    zeros(FT, 51,51)) #0C temperature in all grid cells
atmos = Atmos(grid, 0.0, 0.0, -1.0)  # 0m/s u-velocity, 0.1m/s v-velocity, -1C temperature in all grid cells
# Create domain
domain = Domain( 
  PeriodicBoundary(North, grid), 
  PeriodicBoundary(South, grid), 
  PeriodicBoundary(East, grid),
  PeriodicBoundary(West, grid),
)

coupling_settings = CouplingSettings(
    two_way_coupling_on = true,
    subfloe_point_generator = SubGridPointsGenerator(grid, 1),
)

# Create floes
coord1 = [[[1e-5, 1e5-1e-5], [1e5-1e-5, 1e5-1e-5], [1e5-1e-5, 1e-5], [1e-5, 1e-5], [1e-5, 1e5-1e-5]]]
floe_field = initialize_floe_field(
  FT,
  [coord1],
  domain,
  1.0,  # mean height
  0.0,  # all floes will be the same height
  #rng = Xoshiro(1),
  coupling_settings = coupling_settings,
)

#jldsave("Ocean Stress varfloe3 subgrid/25f_1.0c.jld2"; floe_field)

println(length(floe_field))
# Create model
my_model = Model(grid, ocean, atmos, domain, floe_field) 
# Create simulation
#modulus = 1.5e3*(mean(sqrt.(floe_field.area)) + minimum(sqrt.(floe_field.area))) 
consts = Constants(E = 1.5e8, Cd_ia = 0.0, Cd_ao = 0.0, turnθ = 0.0, f = 0.0)
#println(modulus)
collision_settings = CollisionSettings(
    collisions_on = true,
    floe_floe_max_overlap = 0.55,
    floe_domain_max_overlap = 0.75,
)
fracture_settings = FractureSettings( 
  fractures_on = true,
  criteria = HiblerYieldCurve(floe_field),
  Δt = 75,
  npieces = 3,
  nhistory = 100,
  deform_on = true, 
)


initwriter = InitialStateOutputWriter(
    dir = path,
    filename = "initial_state.jld2",
    overwrite = true,
    )


floewriter = FloeOutputWriter(
    [:coords],
    50,
    dir = path,
    filename = "floes.jld2",
    overwrite = true,
)
gridwriter = GridOutputWriter(
        [:u_grid, :v_grid],  # more options if needed
        50,   # timestep
        grid,
        (51, 51),  # dimension you want the output in... I would just put the dimension of the grid
        dir = path,
        filename = "grid.nc",
        overwrite = true,
        )
writers = OutputWriters(initwriter, floewriter, gridwriter)
simulation = Simulation( 
  model = my_model, 
  consts = consts, 
  Δt = 20,  
  nΔt = 1000, 
  verbose = true, 
  fracture_settings = fracture_settings,
  collision_settings = collision_settings,
  coupling_settings = coupling_settings;
  writers = writers,
)
# Run simulation
function run_with_saving(simulation, tstep_between_save)
    nsave = fld(simulation.nΔt, tstep_between_save) + 1
    τx = zeros(FT, grid.Nx + 1, grid.Ny + 1, nsave)
    uvel = zeros(FT, grid.Nx + 1, grid.Ny + 1, nsave)
    Subzero.startup_sim(simulation, nothing, 1)
    tstep = 0
    counter = 1
    while tstep <= simulation.nΔt
        # Timestep the simulation forward
        timestep_sim!(simulation, tstep)
        if mod(tstep, tstep_between_save) == 0
            τx[:, :, counter] .= simulation.model.ocean.τx
            counter += 1
        end
        tstep+=1
    end
    Subzero.teardown_sim(simulation)
    return τx
end

τx = run_with_saving(simulation, 50)
jldsave("output/coupling_tests/taux.jld2"; τx)