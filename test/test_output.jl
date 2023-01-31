
function test_basic_outputwriters()
    grid = RegRectilinearGrid(-1e5, 1e5, -1e5, 1e5, 1e4, 1e4)
    ocean = Ocean(grid, 0.0, 0.0, 0.0)
    atmos = Atmos(grid, 0.0, 0.0, 0.0)
    domain = Domain(OpenBoundary(grid, North()), OpenBoundary(grid, South()), OpenBoundary(grid, East()), OpenBoundary(grid, West()))
    floe_coords = [[[7.5e4, 7.5e4], [7.5e4, 9.5e4], [9.5e4, 9.5e4], 
                    [9.5e4, 7.5e4], [7.5e4, 7.5e4]]]
    model = Model(grid, ocean, atmos, domain, StructArray([Floe(floe_coords, 0.5, 0.0)]))

    dir = "output/sim"
    simulation = Simulation(model = model, consts = Constants(), nΔt = 500)
    initwriter = InitialStateOutputWriter(dir = dir, filename = "initial_state.jld2", overwrite = true)
    gridwriter = GridOutputWriter(100, grid, (5, 10), dir = dir, filename = "grid.nc", overwrite = true)
    floewriter = FloeOutputWriter([:alive, :coords, :area, :mass, :u, :v], 50, dir = dir, filename = "floe.jld2", overwrite = true)
    checkpointwriter = CheckpointOutputWriter(250, dir = dir, overwrite = true)

    run!(simulation, [initwriter, floewriter, checkpointwriter, gridwriter])

    # Test Initial State Output Writer
    fn = joinpath(dir, "initial_state.jld2")
    file = jldopen(fn, "r")
    @test keys(file) == ["sim"]
    @test typeof(file["sim"]) <: Subzero.Simulation
    close(file)
    rm("output/sim/initial_state.jld2")

    # Test Checkpoint Output Writer
    fn = joinpath(dir, "checkpoint.jld2")
    file = jldopen(fn, "r")
    @test Set(keys(file)) == Set(["ocean", "metadata", "atmos", "floes"])
    @test typeof(file["ocean/0"]) <: Subzero.Ocean
    @test typeof(file["atmos/0"]) <: Subzero.Atmos
    @test eltype(file["floes/0"]) <: Subzero.Floe
    @test keys(file["ocean"]) == ["0", "250", "500"]
    close(file)
    rm(fn)

    # Test Floe Output Writer
    fn = joinpath(dir, "floe.jld2")
    file = jldopen(fn, "r")
    @test Set(keys(file)) ==Set(["alive", "coords", "area", "mass", "u", "v", "metadata"])
    @test length(keys(file["alive"])) == 11
    close(file)
    rm(fn)

    # Test Floe Output Writer
    fn = joinpath(dir, "grid.nc")
    file = Dataset(fn, "r")
    @test Set(keys(file)) == Set(["u_grid", "v_grid", "dudt_grid", "dvdt_grid", "overarea_grid", "mass_grid", "area_grid", "height_grid", "si_frac_grid", "time", "x", "y"])
    @test length(file["x"][:]) == 10
    @test length(file["y"][:]) == 5
    @test length(file["time"][:]) == 6
    close(file)
    rm(fn)
end

@testset "Outputwriters" begin
    test_basic_outputwriters()
end
