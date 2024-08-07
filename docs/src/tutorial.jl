# Let’s run a basic simulation with initially stationary floes pushed into a collision boundary by a uniform, zonally flowing ocean. In this simulation, collisions between floes are on by default and we will enable floe fracturing.  

# Create environment
grid = RegRectilinearGrid(
  (0, 1e5), # xbounds
  (0, 1e5), # ybounds
  1e4,      # grid cell width
  1e4,      # grid cell height
 ) 
ocean = Ocean(grid, 0.25, 0.0, 0.0) # 0.25m/s u-velocity, 0m/s v-velocity, 0C temperature in all grid cells
atmos = Atmos(grid, 0.0, 0.1, -1.0)  # 0m/s u-velocity, 0.1m/s v-velocity, -1C temperature in all grid cells
# Create domain
domain = Domain( 
  CollisionBoundary(North, grid), 
  CollisionBoundary(South, grid), 
  CollisionBoundary(East, grid),
  CollisionBoundary(West, grid),
)
# Create floes
floe_settings = FloeSettings(stress_calculator = DecayAreaScaledCalculator(), min_floe_area = 1e5)
floe_arr = initialize_floe_field(
  Float64,
  100,  # number of floes
  [0.7],  # floe concentration
  domain,
  0.5,  # average floe height
  0.05;  # floe height variability
  floe_settings = floe_settings,
) 
# Create model
model = Model(grid, ocean, atmos, domain, floe_arr) 

# Create settings
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area))) 
consts = Constants(E = modulus) 

fracture_settings = FractureSettings( 
  fractures_on = true,
  criteria = HiblerYieldCurve(floe_arr),
  Δt = 75,
  npieces = 3,
  deform_on = true, 
)
# Create output writers
dir = "output/sim"
initwriter = InitialStateOutputWriter(
    dir = dir,
    filename = "initial_state.jld2",
    overwrite = true,
    )
floewriter = FloeOutputWriter(
    100,
    dir = dir,
    filename = "floes.jld2",
    overwrite = true,
)
writers = OutputWriters(initwriter, floewriter)
# Create simulation
Δt = 10
simulation = Simulation( 
  model = model, 
  consts = consts, 
  Δt = Δt,     # 10 second timesteps 
  nΔt = 10000,  # Run for 10000 timesteps
  verbose = true; 
  fracture_settings = fracture_settings,
  writers = writers,
)
# Run simulation
run!(simulation)

# Plot simulation
plot_sim(
    "output/sim/floes.jld2",
    "output/sim/initial_state.jld2",
    Δt,
    "output/sim/example.mp4",
)