@testset "Bin floes" begin
    # Check assert with zero bins

    # Check only one list with Nx = Ny = 1

    # Check one floe centroid in each bin

    # Check floe centroid on an edge

end

@testset "Weld floes" begin
    grid = RegRectilinearGrid(
        (0, 1e5),
        (0, 1e5),
        1e4,
        1e4,
    )
    domain = Subzero.Domain(
        OpenBoundary(North, grid),
        OpenBoundary(South, grid),
        OpenBoundary(East, grid),
        OpenBoundary(West, grid),
    )
    consts = Constants()
    coupling_settings = CouplingSettings()
    coords = [
        [[[0.0, 0.0], [0.0, 5e4], [6e4, 5e4], [6e4, 0.0], [0.0, 0.0]]],
        [[[4e4, 0.0], [4e4, 5e4], [1e5, 5e4], [1e5, 0.0], [0.0, 0.0]]],
        [[[2e4, 4e4], [2e4, 8e4], [3e4, 8e4], [3e4, 4e4], [2e4, 4e4]]]
    ]
    floe_base = initialize_floe_field(
        Float64,
        coords,
        domain,
        1.0,
        0.0,
    )
    a1, a2, a3 = floe_base.area
    h1, h2, h3 = floe_base.height
    weld_settings = WeldSettings(
        weld_on = true,
        Δts = [100, 250, 700],
        Nxs = [1, 2, 1],
        Nys = [1, 2, 2],
        max_weld_area = 1e10,  # entire domain
        welding_coeff = 1000  # will weld for sure
    )
    # Test no floes welding because centroids in different bins (Nx = 2, Ny = 2)
    floes = deepcopy(floe_base)
    Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        coupling_settings,
        weld_settings,
        2, # second set in weld setting lists (Nx = 2, Ny = 2)
        consts,
        10,  # Δt
    )
    @test Subzero.active == floes.status[1].tag ==
        floes.status[2].tag == floes.status[3].tag
    @test floes.area[1] == a1 && floes.area[2] == a2 && floes.area[3] == a3
    @test floes.height[1] == h1 && floes.height[2] == h2 &&
        floes.height[3] == h3

    # Test two floes welding with centroids in the same bin (Nx = 1, Ny = 2)
    floes = deepcopy(floe_base)
    Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        coupling_settings,
        weld_settings,
        3, # third set in weld setting lists (Nx = 1, Ny = 2)
        consts,
        10,  # Δt
    )
    @test Subzero.active == floes.status[1].tag == floes.status[3].tag
    @test Subzero.remove == floes.status[2].tag
    @test floes.area[1] > a1
    @test floes.area[1] == 5e9
    @test floes.area[3] == a3
    @test floes.height[1] > h1
    @test floes.height[3] == h3

    # Test all three floes welding together (Nx = 1, Ny = 1)
    floes = deepcopy(floe_base)
    Subzero.timestep_welding!(
        floes,
        maximum(floes.id),
        grid,
        coupling_settings,
        weld_settings,
        1, # third set in weld setting lists (Nx = 1, Ny = 1)
        consts,
        10,  # Δt
    )
    @test Subzero.active == floes.status[1].tag
    @test Subzero.remove == floes.status[2].tag == floes.status[3].tag
    @test floes.area[1] == 5.3e9
    @test floes.area[1] > a1
    @test floes.height[1] > h1

    # All floes are too big to weld

    # All floes are too small to weld

    # Floe only welds with largest interaction (small max weld area)

end