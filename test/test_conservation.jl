function conservation_simulation(grid, domain, floes)
    ocean = Ocean(grid, 0.0, 0.0, 0.0)
    atmos = Atmos(grid, 0.0, 0.0, 0.0)
    model = Model(grid, ocean, atmos, domain, floes)
    dir = "output/conservation"
    initwriter = InitialStateOutputWriter(
        dir = dir,
        overwrite = true,
    )
    floewriter = FloeOutputWriter(
        10,
        dir = dir,
        overwrite = true,
    )
    writers = OutputWriters(
        initialwriters = StructArray([initwriter]),
        floewriters = StructArray([floewriter]),
    )
    modulus = 1.5e3*(mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area)))
    # No friction
    consts = Constants(E = modulus, μ = 0.0)
    # No ocean/atmosphere interactions
    coupling_settings = CouplingSettings(coupling_on = false)

    simulation = Simulation(
        model = model,
        consts = consts,
        Δt = 10,
        nΔt = 10000,
        verbose = false,
        coupling_settings = coupling_settings,
        writers = writers,
    )
    run!(simulation)
    em_lists = check_energy_momentum_conservation_julia(
        joinpath(dir, "floes.jld2"),
        dir,
    )
    # Note that if any of these start with 0, this won't work
    return [vals[1] != 0 ?
        (vals[end] - vals[1])/vals[1] * 100 :
        @warn "Starting value is 0 for series, can't calculate percent change"
        for vals in em_lists]
end


@testset "Conservation of Energy and Momentum" begin
    grid = RegRectilinearGrid(
        Float64,
        (0, 1e5),
        (0, 1e5),
        1e4,
        1e4,
    )
    collision_domain = Domain(CollisionBoundary(grid, North()),
        CollisionBoundary(grid, South()),
        CollisionBoundary(grid, East()),
        CollisionBoundary(grid, West()),
    )
    open_domain = Domain(OpenBoundary(grid, North()),
        OpenBoundary(grid, South()),
        OpenBoundary(grid, East()),
        OpenBoundary(grid, West()),
    )
    rng = Xoshiro(1)
    floe1 = [[[2e4, 2e4], [2e4, 5e4], [5e4, 5e4], [5e4, 2e4], [2e4, 2e4]]]
    floe2 = [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]
    floe3 = [[[5.5e4, 2e4], [5.25e4, 4e4], [5.75e4, 4e4], [5.5e4, 2e4]]]
    # One block hits the wall
    floe_on_wall = initialize_floe_field(
        [Subzero.translate(floe1, [-1.9e4, 0.0])],
        open_domain, # Just affects shape, type doesnt' matter
        0.25,
        0.0,
        rng = rng,
        nhistory = 100,
    )
    floe_on_wall.u[1] = -0.1
    floe_on_wall.v[1] = 0.01
    @test all(abs.(conservation_simulation(grid, open_domain, floe_on_wall)) .< 1)
    # Two blocks crashing head on - no rotation
    head_on_floes = initialize_floe_field(
        [floe1, floe2],
        open_domain, # Just affects shape, type doesnt' matter
        0.25,
        0.0,
        rng = rng,
        nhistory = 100,
    )
    head_on_floes.u[1] = 0.15
    head_on_floes.u[2] = -0.1
    head_on_floes.v[1] = 0.02
    head_on_floes.v[2] = 0.02
    head_on_floes.ξ[1] = 1e-7
    # less than 1% change in both with open domain
    @test all(abs.(conservation_simulation(grid, open_domain, head_on_floes)) .< 1)

    # Two blocks crashing offset - rotation
    offset_floes = initialize_floe_field(
        [floe1, Subzero.translate(floe2, [0.0, 1e4])],
        open_domain, # Just affects shape, type doesnt' matter
        0.25,
        0.0,
        rng = rng,
        nhistory = 100,
    )
    offset_floes.u[1] = 0.11
    offset_floes.u[2] = -0.1
    offset_floes.v[1] = 0.02
    offset_floes.v[2] = 0.02
    offset_floes.ξ[1] = 1e-7
    @test all(abs.(conservation_simulation(grid, open_domain, offset_floes)) .< 1)

    # Two rectangular boxes with a triangle inbetween causing rotation
    rotating_floes = initialize_floe_field(
        [floe1, floe2, floe3],
        open_domain, # Just affects shape, type doesnt' matter
        0.25,
        0.0,
        rng = rng,
        nhistory = 100,
    )
    rotating_floes.u[1] = 0.11
    rotating_floes.u[2] = -0.1
    rotating_floes.ξ[3] = 1e-5
    rotating_floes.v[1] = 0.001
    rotating_floes.v[2] = 0.001
    rotating_floes.v[3] = 0.001
    @test all(abs.(conservation_simulation(grid, open_domain, rotating_floes)) .< 1)
end