using LibGEOS
@testset "Ridging and Rafting" begin
    function update_height(floes, i, new_height, consts)
        floes.height[i] = new_height
        floes.mass[i] = floes.area[i] * floes.height[i] * consts.ρi
        floes.moment[i] = Subzero.calc_moment_inertia(
            floes.coords[i],
            floes.centroid[i],
            floes.height[i],
            ρi = consts.ρi,
        )
    end

    function test_ridge_raft_scenario(
        rr_settings,
        new_height1,
        new_height2,
        floe1_subsume,
        floes,
        domain,
        coupling_settings,
        simp_settings,
        consts,
    )
        # Test floe 1 > min ridge height and floe 2 < min ridge height
        update_height(floes, 1, new_height1, consts)
        update_height(floes, 2, new_height2, consts)
        mass1, mass2 = floes.mass
        moment1, moment2 = floes.moment
        height1, height2 = floes.height
        cent1, cent2 = floes.centroid
        Subzero.timestep_ridging_rafting!(
            floes,
            2,
            domain,
            rr_settings,
            coupling_settings,
            simp_settings,
            consts,
            10,
        )
        @test mass1 + mass2 == sum(floes.mass)
        @test floe1_subsume ?
            (mass1 < floes.mass[1] && mass2 > floes.mass[2]) :
            (mass1 > floes.mass[1] && mass2 < floes.mass[2])
        @test floe1_subsume ?
            (floes.height[1] > height1 && floes.height[2] == height2) :
            (floes.height[1] == height1 && floes.height[2] > height2)
        @test floe1_subsume ?
            moment1 * floes.height[1] / height1 == floes.moment[1] :
            moment2 * floes.height[2] / height2 == floes.moment[2]
        @test floe1_subsume ?
            cent1 == floes.centroid[1] && cent2 != floes.centroid[2] :
            cent2 == floes.centroid[1] && cent1 != floes.centroid[1]
        @test LibGEOS.area(LibGEOS.intersection(
            LibGEOS.Polygon(floes.coords[1]),
            LibGEOS.Polygon(floes.coords[2])
        )) == 0  # floes DO NOT overlap anymore!
    end

    # Setup for tests
    grid = RegRectilinearGrid(
           (0, 1e5),
           (0, 1e5),
           1e4,
           1e4,
    )
    domain = Subzero.Domain(
           CollisionBoundary(North, grid),
           CollisionBoundary(South, grid),
           CollisionBoundary(East, grid),
           CollisionBoundary(West, grid),
    )
    consts = Constants()
    coupling_settings = CouplingSettings()
    simp_settings = SimplificationSettings()
    collision_settings = CollisionSettings(floe_floe_max_overlap = 0.99) # don't fuse
    lock = Threads.SpinLock()
    coords = [
        [[[0.1e4, 0.1e4], [0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [0.1e4, 0.1e4]]],
        [[[1.8e4, 1.8e4], [1.8e4, 4e4], [4e4, 4e4], [4e4, 1.8e4], [1.8e4, 1.8e4]]],
    ]
    floes = initialize_floe_field(
        Float64,
        coords,
        domain,
        1.0,
        0.0,
    )
    Subzero.timestep_collisions!(  # Add interactions
        floes,
        2,
        domain,
        consts,
        10,
        collision_settings, 
        lock,
    )
    @testset "Ridge Floe-Floe" begin
        rr_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # force ridging
            raft_probability = 0.0,
        )
        # Test floe 2 ridging on top of floe 1
        test_ridge_raft_scenario(
            rr_settings,
            1.0,
            0.2,
            true,  # floe 1 should subsume floe 2
            floes,
            domain,
            coupling_settings,
            simp_settings,
            consts,
        )
        # Test floe 1 ridging on top of floe 2
        # test_ridge_raft_scenario(
        #     rr_settings,
        #     0.2,
        #     0.1,
        #     false,  # floe 2 should subsume floe 1
        #     floes,
        #     domain,
        #     coupling_settings,
        #     simp_settings,
        #     consts,
        # )
    end
    @testset "Ridge Floe-Domain" begin
        
    end
    @testset "Raft Floe-Floe" begin
        rr_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 0.0,  # force rafting
            raft_probability = 1.0,
        )
        
    end
    @testset "Raft Floe-Domain" begin
        
    end 
end