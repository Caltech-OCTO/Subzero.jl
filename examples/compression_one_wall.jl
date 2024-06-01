using JLD2, Random, Statistics, Subzero, DelimitedFiles

# User Inputs
const FT = Float64
const Lx = 1e5
const Ly = 1e5
const Δgrid = 2e3
const hmean = 0.25
const Δh = 0.125
const Δt = 20 # seconds
# const theta = 2*pi/3

# Model instantiation
grid = RegRectilinearGrid(
    (0.0, Lx),
    (0.0, Ly),
    Δgrid,
    Δgrid,
)
ocean = Ocean(grid, 0.0, 0.0, 0.0)
atmos = Atmos(grid, 0.0, 0.0, -1.0)


# V = 0.2
dir = "output/moving_wall"
# for theta in [pi/2, pi/3, pi/4, pi/6]
for theta in [pi/2]
    Vs = 0.005:0.005:0.1
    tbs = []
    # for V in Vs
    for V in [0.01]
        # Domain creation
        nboundary = CollisionBoundary(North, grid)
        sboundary = MovingBoundary(South, grid, 0.0, V*sin(theta))
        eboundary = CollisionBoundary(East, grid)
        wboundary = MovingBoundary(West, grid, V*cos(theta), 0.0)

        domain = Domain(nboundary, sboundary, eboundary, wboundary)

        # Floe creation
        floe_arr = initialize_floe_field(
            FT,
            # [[[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]],
            [[[[0, 0], [Lx, 0], [Lx, Ly], [0, Ly], [0, 0]]]],
            domain,
            hmean,  # mean height of 0.25
            0.0,    # all floes will be the same height
            Δt;  
            rng = Xoshiro(1),
        )
        # display(floe_arr.area)
        nfloes = length(floe_arr)
        floe_arr.u .= 0
        floe_arr.v .= -0.01
        # Model creation
        model = Model(grid, ocean, atmos, domain, floe_arr)

        # Simulation setup
        modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
        consts = Constants(E = modulus, Cd_io = 0.0, f = 0.0, turnθ = 0.0, Cd_ia = 0.0)

        # Run simulation
        # run_until_break!
        run_time!(simulation) = @time run!(simulation)
        # run_time_till_break!(simulation) = @time run_until_break!(simulation)
        

        # Output setup
        initwriter = InitialStateOutputWriter(dir = dir, overwrite = true)
        floewriter = FloeOutputWriter(40, dir = dir, overwrite = true)
        writers = OutputWriters(initwriter, floewriter)

        # Simulation settings 

        coupling_settings = CouplingSettings(two_way_coupling_on = false)

        fracture_settings = FractureSettings(
                fractures_on = true,
                criteria = HiblerYieldCurve(floe_arr),
                Δt = 1, # in timesteps
                npieces = 3,
                deform_on = false,
        )

        simulation = Simulation(
            model = model,
            consts = consts,
            Δt = Δt,
            nΔt = 1000,
            verbose = true,
            writers = writers,
            rng = Xoshiro(1),
            coupling_settings = coupling_settings, 
            fracture_settings = fracture_settings
        )


        steps_till_break = run_time!(simulation)
        time_till_break = steps_till_break*Δt

        # run_time_till_break!(simulation)

        # function run_until_break!(sim; logger = nothing, messages_per_tstep = 1)
        #     startup_sim(sim, logger, messages_per_tstep)
        #     tstep = 0
        #     while length(floe_arr <= 1)
        #         # Timestep the simulation forward
        #         timestep_sim!(sim, tstep)
        #         tstep+=1
        #     end
        #     teardown_sim(sim)
        #     return
        # end

        # 

        Subzero.plot_sim(
            dir*"/floes.jld2",
            dir*"/initial_state.jld2",
            20,
            dir*"/TEST.mp4",
        )
        push!(tbs, time_till_break)
    end

    # writedlm( dir*"/tbs_theta_$(theta).csv", tbs, ',')
end


## multiply tstep by time step to get in seconds
## can make plots of different values of theta and V
## can populate array with all the timesteps, and then write to a csv or jld2