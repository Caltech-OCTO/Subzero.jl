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

    function assign_random_velocities!(floes)
        for i in eachindex(floes)
            floes.u[i] = (-1)^rand(0:1) * rand()
            floes.v[i] = (-1)^rand(0:1) * rand()
            floes.ξ[i] = (-1)^rand(0:1) * 0.05rand()
            floes.p_dxdt[i] = (-1)^rand(0:1) * rand()
            floes.p_dydt[i] = (-1)^rand(0:1) * rand()
            floes.p_dαdt[i] = (-1)^rand(0:1) * 0.05rand()
        end
        return
    end
    function setup_floes_with_inters(coords, domain, consts,
        collision_settings, lock,  Δx = nothing, Δy = nothing,
    )
        floes = initialize_floe_field(
            Float64,
            coords,
            domain,
            1.0,
            0.0,
        )
        if !isnothing(Δx)
            for i in eachindex(Δx)
                Subzero.translate!(floes.coords[i], Δx[i], Δy[i])
                floes.centroid[i][1] += Δx[i]
                floes.centroid[i][2] += Δy[i]
            end
        end
        assign_random_velocities!(floes)
        add_ghosts!(floes, domain)
        Subzero.timestep_collisions!(  # Add interactions
            floes,
            length(floes),
            domain,
            consts,
            10,
            collision_settings, 
            lock,
        )
        return floes
    end
    function calc_needed_momentum(floes)
        x_momentum, y_momentum = Subzero.calc_linear_momentum(
            floes.u,
            floes.v,
            floes.mass,
        )
        spin_momentum, angular_momentum = Subzero.calc_angular_momentum(
            floes.u,
            floes.v,
            floes.mass,
            floes.ξ,
            floes.moment,
            first.(floes.centroid),
            last.(floes.centroid),
        )
        p_x_momentum, p_y_momentum = Subzero.calc_linear_momentum(
            floes.p_dxdt,
            floes.p_dydt,
            floes.mass,
        )
        p_spin_momentum, p_orbital_momentum = Subzero.calc_angular_momentum(
            floes.p_dxdt,
            floes.p_dydt,
            floes.mass,
            floes.p_dαdt,
            floes.moment,
            first.(floes.centroid) .* floes.p_dxdt,
            last.(floes.centroid) .* floes.p_dydt,
        )
        return x_momentum, y_momentum, spin_momentum, angular_momentum,
            p_x_momentum, p_y_momentum, p_spin_momentum, p_orbital_momentum
    end

    function conservation_of_momentum_tests(floes,
        x_momentum_init, y_momentum_init,
        spin_momentum_init, orbital_momentum_init,
        p_x_momentum_init, p_y_momentum_init,
        p_spin_momentum_init, p_orbital_momentum_init,
    )
        x_momentum_after, y_momentum_after,
        spin_momentum_after, orbital_momentum_after,
        p_x_momentum_after, p_y_momentum_after,
        p_spin_momentum_after, p_orbital_momentum_after = calc_needed_momentum(floes)

        @test isapprox(x_momentum_init, x_momentum_after, atol = 1e-3)
        @test isapprox(y_momentum_init, y_momentum_after, atol = 1e-3)
        @test isapprox(p_x_momentum_init, p_x_momentum_after, atol = 1e-3)
        @test isapprox(p_y_momentum_init, p_y_momentum_after, atol = 1e-3)
        @test isapprox(
            spin_momentum_init + orbital_momentum_init,
            spin_momentum_after + orbital_momentum_after,
            atol = 1e3,
        )
        @test isapprox(
            p_spin_momentum_init + p_orbital_momentum_init,
            p_spin_momentum_after + p_orbital_momentum_after,
            atol = 1e3,
        )
    end

    function test_floe_floe_rr_scenario(
        rr_settings,
        new_height1,
        new_height2,
        floe1_subsume,
        floe2_subsume,
        floes,
        domain,
        coupling_settings,
        simp_settings,
        consts,
    )
        @assert !floe1_subsume || !floe2_subsume "Both floe 1 and floe 2 can't subsume for testing"
        # Test floe 1 > min ridge height and floe 2 < min ridge height
        update_height(floes, 1, new_height1, consts)
        update_height(floes, 2, new_height2, consts)
        mass1, mass2 = floes.mass
        moment1, moment2 = floes.moment
        height1, height2 = floes.height
        cent1, cent2 = floes.centroid
        # initial momentum
        x_momentum_init, y_momentum_init,
        spin_momentum_init, orbital_momentum_init,
        p_x_momentum_init, p_y_momentum_init,
        p_spin_momentum_init, p_orbital_momentum_init = calc_needed_momentum(floes)
        # Ridge and raft floes
        Subzero.timestep_ridging_rafting!(
            floes,
            StructArray{Floe{Float64}}(undef, 0),
            2,
            domain,
            maximum(floes.id),
            rr_settings,
            coupling_settings,
            simp_settings,
            consts,
            10,
        )
        @test mass1 + mass2 == sum(floes.mass)
        conservation_of_momentum_tests(floes,
            x_momentum_init, y_momentum_init,
            spin_momentum_init, orbital_momentum_init,
            p_x_momentum_init, p_y_momentum_init,
            p_spin_momentum_init, p_orbital_momentum_init,
        )
        if floe1_subsume || floe2_subsume
            @test LibGEOS.area(LibGEOS.intersection(
                LibGEOS.Polygon(floes.coords[1]),
                LibGEOS.Polygon(floes.coords[2])
            )) == 0  # floes DO NOT overlap anymore!
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
                cent2 == floes.centroid[2] && cent1 != floes.centroid[1]
        else
            @test mass1 == floes.mass[1] && mass2 == floes.mass[2]
            @test floes.height[1] == height1 && floes.height[2] == height2
            @test moment1 == floes.moment[1] && moment2 == floes.moment[2]
            @test cent1 == floes.centroid[1] && cent2 == floes.centroid[2]
        end
    end
    function test_floe_domain_rr_scenario(rr_settings, does_raft, floes,
        domain, boundary_poly, topo_poly, bounds_overlap_area, topo_overlap_area,
        height1, height2, coupling_settings, simp_settings, consts)
        update_height(floes, 1, height1, consts)
        update_height(floes, 2, height2, consts)
        total_mass = sum(floes.mass)
        h1, h2 = floes.height
        area1, area2 = floes.area
        cent1, cent2 = floes.centroid

        x_momentum_init, y_momentum_init,
        spin_momentum_init, orbital_momentum_init,
        p_x_momentum_init, p_y_momentum_init,
        p_spin_momentum_init, p_orbital_momentum_init = calc_needed_momentum(floes)
        Subzero.timestep_ridging_rafting!(
            floes,
            StructArray{Floe{Float64}}(undef, 0),
            2,
            domain,
            maximum(floes.id),
            rr_settings,
            coupling_settings,
            simp_settings,
            consts,
            10,
        )
        conservation_of_momentum_tests(floes,
            x_momentum_init, y_momentum_init,
            spin_momentum_init, orbital_momentum_init,
            p_x_momentum_init, p_y_momentum_init,
            p_spin_momentum_init, p_orbital_momentum_init,
        )
        if does_raft
            @test total_mass > sum(floes.mass)
            @test h1 == floes.height[1] && h2 == floes.height[2]
            @test area1 - bounds_overlap_area == floes.area[1]
            @test area2 - topo_overlap_area == floes.area[2]
            @test cent1 != floes.centroid[1] && cent2 != floes.centroid[2]
            @test LibGEOS.area(LibGEOS.intersection(
                    LibGEOS.Polygon(floes.coords[1]),
                    boundary_poly,
            )) == 0
            @test LibGEOS.area(LibGEOS.intersection(
                LibGEOS.Polygon(floes.coords[2]),
                topo_poly,
            )) == 0
        else
            @test total_mass == sum(floes.mass)
            @test h1 == floes.height[1] && h2 == floes.height[2]
            @test area1 == floes.area[1]
            @test area2 == floes.area[2]
            @test cent1 == floes.centroid[1] && cent2 == floes.centroid[2]
        end
    end

    function test_floe_ghost_rr_scenario(rr_settings, floes, height1, height2,
        floe1_subsume, domain, coupling_settings, simp_settings, consts,
    )
        update_height(floes, 1, height1, consts)  # floe 1 will ridge onto floe 2
        update_height(floes, 3, height1, consts)  # floe 1 has ghost floe 3
        update_height(floes, 2, height2, consts)
        total_mass = floes.mass[1] + floes.mass[2]
        h1, h2, h3 = floes.height
        area1, area2, area3 = floes.area
        cent1, cent2, cent3 = floes.centroid
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            StructArray{Floe{Float64}}(undef, 0),
            2,
            domain,
            maximum(floes.id),
            rr_settings,
            coupling_settings,
            simp_settings,
            consts,
            10,
        )
        @test length(floes) == 3
        # conservation of mass doesn't include ghost floes bc they're not real
        @test total_mass == floes.mass[1] + floes.mass[2]
        @test floe1_subsume ?
            h1 < floes.height[1] && h3 < floes.height[3] &&  h1 == h3 :
            h1 == floes.height[1] && h3 == floes.height[3] && h1 == h3
        @test floe1_subsume ? 
            h2 == floes.height[2] :
            h2 < floes.height[2]
        @test floe1_subsume ?
            area1 == floes.area[1] && area3 == floes.area[3] && area1 == area3 :
            area1 > floes.area[1] && area3 > floes.area[3] && area1 == area3
        @test floe1_subsume ?
            area2 > floes.area[2] :
            area2 == floes.area[2]
        @test floe1_subsume ?
            cent2 != floes.centroid[2] :
            cent2 == floes.centroid[2]
        @test floe1_subsume ?
            cent1 == floes.centroid[1] && cent3 == floes.centroid[3] :
            cent1 != floes.centroid[1] && cent3 != floes.centroid[3]
        @test floes.centroid[1] .- floes.centroid[3] == cent1 .- cent3
    end

    # Setup for tests
    grid = RegRectilinearGrid(
           (0, 1e5),
           (0, 1e5),
           1e4,
           1e4,
    )
    topo_coords = [[[5e4, 5e4], [5e4, 7e4], [7e4, 7e4], [7e4, 5e4], [5e4, 5e4]]]
    collision_domain = Subzero.Domain(
        CollisionBoundary(North, grid),
        CollisionBoundary(South, grid),
        CollisionBoundary(East, grid),
        CollisionBoundary(West, grid),
        initialize_topography_field([topo_coords])
    )
    periodic_domain = Subzero.Domain(
        PeriodicBoundary(North, grid),
        PeriodicBoundary(South, grid),
        PeriodicBoundary(East, grid),
        PeriodicBoundary(West, grid),
    )
    boundary_poly = LibGEOS.union(
        LibGEOS.MultiPolygon([collision_domain.north.coords, collision_domain.south.coords]),
        LibGEOS.MultiPolygon([collision_domain.east.coords, collision_domain.west.coords]),
    )
    topo_poly = LibGEOS.Polygon(topo_coords)
    consts = Constants()
    coupling_settings = CouplingSettings()
    simp_settings = SimplificationSettings()
    collision_settings = CollisionSettings(floe_floe_max_overlap = 0.99) # don't fuse
    lock = Threads.SpinLock()
    @testset "Floe-Floe Ridging and Rafting" begin
        coords = [
            [[[0.1e4, 0.1e4], [0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [0.1e4, 0.1e4]]],
            [[[1.8e4, 1.8e4], [1.8e4, 4e4], [4e4, 4e4], [4e4, 1.8e4], [1.8e4, 1.8e4]]],
        ]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock,
        )
        # Test scenario with no ridging or rafting
        no_rr_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 0.0,  # no ridging
            raft_probability = 0.0,  # no rafting
        )
        floes = deepcopy(floes_base)
        # Test floe 2 ridging on top of floe 1
        test_floe_floe_rr_scenario(
            no_rr_settings,
            1.0,
            1.0,
            false,
            false,
            floes,
            collision_domain,
            coupling_settings,
            simp_settings,
            consts,
        )
        # Ridging Tests
        ridge_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # force ridging
            raft_probability = 0.0,
        )
        floes = deepcopy(floes_base)
        # Test floe 2 ridging on top of floe 1
        println("2")
        test_floe_floe_rr_scenario(
            ridge_settings,
            1.0,
            0.1,
            true,  # floe 1 should subsume floe 2
            false,
            floes,
            collision_domain,
            coupling_settings,
            simp_settings,
            consts,
        )
        floes = deepcopy(floes_base)
        # Test floe 1 ridging on top of floe 2
        println("3")
        test_floe_floe_rr_scenario(
            ridge_settings,
            0.1,
            1.0,
            false,
            true,
            floes,
            collision_domain,
            coupling_settings,
            simp_settings,
            consts,
        )
        floes = deepcopy(floes_base)
        # Test both floes being too thin to ridge
        println("4")
        test_floe_floe_rr_scenario(
            ridge_settings,
            0.1,
            0.1,
            false,
            false,
            floes,
            collision_domain,
            coupling_settings,
            simp_settings,
            consts,
        )
        # Rafting Tests
        raft_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 0.0,  # force rafting
            raft_probability = 1.0,
            max_floe_raft_height = 1.0,
        )
        floes = deepcopy(floes_base)
        # Test floe 2 rafting on top of floe 1
        println("5")
        test_floe_floe_rr_scenario(
            raft_settings,
            1.0,
            0.001,  # bias so that floe 1 will subsume floe 2
            true,  # floe 1 should subsume floe 2
            false,
            floes,
            collision_domain,
            coupling_settings,
            simp_settings,
            consts,
        )
        # Test floe 2 rafting on top of floe 1
        floes = deepcopy(floes_base)
        println("6")
        test_floe_floe_rr_scenario(
            raft_settings,
            0.001, # bias so that floe 2 will subsume floe 1
            1.0,
            false,
            true,  # floe 2 should subsume floe 1
            floes,
            collision_domain,
            coupling_settings,
            simp_settings,
            consts,
        )
    end
    @testset "Floe-Domain Ridging and Rafting" begin
        coords = [
            [[[-0.1e4, -0.1e4], [-0.1e4, 2e4], [2e4, 2e4], [2e4, -0.1e4], [-0.1e4, -0.1e4]]],
            [[[3e4, 3e4], [3e4, 5e4], [5e4, 5e4], [5e4, 3e4], [3e4, 3e4]]]
        ]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock, [0.0, 0.5e4], [0.0, 0.5e4],
        )
        bounds_overlap_area = LibGEOS.area(LibGEOS.intersection(
                LibGEOS.Polygon(floes_base.coords[1]),
                boundary_poly,
        ))
        topo_overlap_area = LibGEOS.area(LibGEOS.intersection(
                LibGEOS.Polygon(floes_base.coords[2]),
                topo_poly,
        ))
        # Ridging with domain
        floes = deepcopy(floes_base)
        println("1")
        ridge_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # force ridging
            raft_probability = 0.0,
        )
        test_floe_domain_rr_scenario(ridge_settings, true, floes, collision_domain,
            boundary_poly, topo_poly, bounds_overlap_area, topo_overlap_area,
            0.1, 0.1, coupling_settings, simp_settings, consts,
        )
        # Not ridging with domain
        println("2")
        floes = deepcopy(floes_base)
        test_floe_domain_rr_scenario(ridge_settings, false, floes, collision_domain,
            boundary_poly, topo_poly, bounds_overlap_area, topo_overlap_area, 
            2.0, 2.0, coupling_settings, simp_settings, consts,
        )
        # Rafting with domain
        floes = deepcopy(floes_base)
        rafting_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 0.0,  # force ridging
            raft_probability = 1.0,
        )
        println("3")
        test_floe_domain_rr_scenario(rafting_settings, true, floes, collision_domain,
            boundary_poly, topo_poly, bounds_overlap_area, topo_overlap_area, 
            0.1, 0.1, coupling_settings, simp_settings, consts,
        )
        # Not rafting with domain   
        floes = deepcopy(floes_base)
        println("4")
        test_floe_domain_rr_scenario(rafting_settings, false, floes, collision_domain,
            boundary_poly, topo_poly, bounds_overlap_area, topo_overlap_area, 
            0.3, 0.3, coupling_settings, simp_settings, consts,
        )
    end

    @testset "Special Ridge Raft Cases" begin
        ridge_settings = Subzero.RidgeRaftSettings(
            ridge_probability = 1.0,  # no ridging
            raft_probability = 0.0,  # no rafting
        )
        # first floe overlaps with both other floes, and it will break when ridged with floe 2
        coords = [
            [[[2.75e4, 0.75e4], [0.75e4, 2.75e4], [1.25e4, 2.75e4], [3.25e4, 0.75e4], [2.75e4, 0.75e4]]],
            [[[0.1e4, 0.1e4], [0.1e4, 2.25e4], [2.25e4, 2.25e4], [2.25e4, 0.1e4], [0.1e4, 0.1e4]]],
            [[[2.5e4, 0.1e4], [2.5e4, 2.25e4], [3e4, 2.25e4], [3e4, 0.1e4], [2.5e4, 0.1e4]]],
        ]
        floes_base = setup_floes_with_inters(coords, collision_domain, consts,
            collision_settings, lock,
        )
        # Test scenario with floe breaking into pieces
        floes = deepcopy(floes_base)
        update_height(floes, 1, 0.1, consts)  # floe1 will ridge onto floe 2
        total_mass = sum(floes.mass)
        h1, h2, h3 = floes.height
        area1, area2, area3 = floes.area
        cent1, cent2, cent3 = floes.centroid
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            3,
            collision_domain,
            maximum(floes.id),
            ridge_settings,
            coupling_settings,
            simp_settings,
            consts,
            10,
        )
        @test length(pieces_list) == 1
        @test total_mass == sum(floes.mass) + sum(pieces_list.mass)
        # Make sure floes 1 ridged onto floe 2 and broke
        @test h1 == floes.height[1] && h2 < floes.height[2]
        @test floes.height[1] == pieces_list.height[1]
        @test cent2 == floes.centroid[2]
        # Make sure floe 3 wasn't ridged/rafted
        @test h3 == floes.height[3]
        @test cent3 == floes.centroid[3]
        # Make sure IDs are correct
        @test max_id == 5
        @test floes.id == [4, 2, 3] && pieces_list.id[1] == 5
        @test floes.parent_ids[1] == [1] && pieces_list.parent_ids[1] == [1]
        @test isempty(floes.parent_ids[2]) && isempty(floes.parent_ids[3])

        # Test parent ridging -> update ghost floes
        coords = [
            [[[-0.1e4, 0.1e4], [-0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [-0.1e4, 0.1e4]]],
            [[[1.8e4, 1.8e4], [1.8e4, 4e4], [4e4, 4e4], [4e4, 1.8e4], [1.8e4, 1.8e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
        # Test parent-parent ridge, no breakage and update ghost
            # parent with ghost is subsumed by floe 2
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 0.1, 1.0,
            false, periodic_domain, coupling_settings, simp_settings, consts,
        )
            # parent with ghost subsumes floe 2 
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 1.0, 0.1,
            true, periodic_domain, coupling_settings, simp_settings, consts,
        )

        # Test parent-ghost ridge, no breakage and update parents
        coords = [
            [[[-0.1e4, 0.1e4], [-0.1e4, 2e4], [2e4, 2e4], [2e4, 0.1e4], [-0.1e4, 0.1e4]]],
            [[[8e4, 1.8e4], [8e4, 4e4], [9.92e4, 4e4], [9.92e4, 1.8e4], [8e4, 1.8e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
            # ghost (floe 3) is subsumed by floe 2
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 0.1, 1.0,
            false, periodic_domain, coupling_settings, simp_settings, consts,
        )
            # floe 2 is subsumed by ghost (floe 3)
        floes = deepcopy(base_floes)
        test_floe_ghost_rr_scenario(ridge_settings, floes, 1.0, 0.1,
            true, periodic_domain, coupling_settings, simp_settings, consts,
        )
        
        # Test parent-parent ridge, breakage
        coords = [
            [[[-1e4, 0.75e4], [0.75e4, 2.75e4], [1.25e4, 2.75e4], [-0.5e4, 0.75e4], [-1e4, 0.75e4]]],
            [[[0.1e4, 0.1e4], [0.1e4, 2.25e4], [2.25e4, 2.25e4], [2.25e4, 0.1e4], [0.1e4, 0.1e4]]],
            [[[9e4, 0.1e4], [9e4, 2.25e4], [9.5e4, 2.25e4], [9.5e4, 0.1e4], [9e4, 0.1e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
        floes = deepcopy(base_floes)
        update_height(floes, 1, 0.1, consts)
        update_height(floes, 4, 0.1, consts)  # ghost of floe 1
        total_mass = sum(floes.mass[1:3])
        h1, h2, h3, h4 = floes.height
        area1, area2, area3, area4 = floes.area
        cent1, cent2, cent3, cent4 = floes.centroid
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            3,
            periodic_domain,
            maximum(floes.id),
            ridge_settings,
            coupling_settings,
            simp_settings,
            consts,
            10,
        )
        # Make sure each piece is saved, ghosts are marked to remove and
        # don't interact with any other floes
        @test length(pieces_list) == 1
        @test total_mass == sum(floes.mass[1:3]) + sum(pieces_list.mass)
        # Make sure floes 1 ridged onto floe 2 and broke
        @test isapprox(h1, floes.height[1], atol = 1e-15) && h2 < floes.height[2]
        @test floes.height[1] == pieces_list.height[1]
        @test cent2 == floes.centroid[2]
        # Make sure floe 3 wasn't ridged/rafted
        @test h3 == floes.height[3]
        @test cent3 == floes.centroid[3]
        # Make sure IDs are correct
        @test max_id == 5
        @test floes.id == [4, 2, 3, 1] && pieces_list.id[1] == 5
        @test floes.parent_ids[1] == [1] && pieces_list.parent_ids[1] == [1]
        @test isempty(floes.parent_ids[2]) && isempty(floes.parent_ids[3])
        # Make sure ghost floe isn't changed and is marked for remove_floe_overlap
        @test floes.status[4].tag == Subzero.remove
        @test h4 == floes.height[4]
        @test cent4 == floes.centroid[4]
        @test area4 == floes.area[4]

        # Test parent-ghost ridge, breakage
            # Make sure each piece is saved, ghosts are marked to remove and
            # don't interact with any other floes
        coords = [
            [[[8.5e4, 0.75e4], [10.75e4, 2.75e4], [11.25e4, 2.75e4], [9e4, 0.75e4], [8.5e4, 0.75e4]]],
            [[[8.5e4, 0.1e4], [8.5e4, 2.25e4], [9e4, 2.25e4], [9e4, 0.1e4], [8.5e4, 0.1e4]]],
            [[[0.1e4, 0.1e4], [0.1e4, 2.25e4], [2.25e4, 2.25e4], [2.25e4, 0.1e4], [0.1e4, 0.1e4]]],
        ]
        base_floes = setup_floes_with_inters(coords, periodic_domain, consts,
            collision_settings, lock
        )
        floes = deepcopy(base_floes)
        update_height(floes, 1, 0.1, consts)
        update_height(floes, 4, 0.1, consts)  # ghost of floe 1
        total_mass = sum(floes.mass[1:3])
        h1, h2, h3, h4 = floes.height
        area1, area2, area3, area4 = floes.area
        cent1, cent2, cent3, cent4 = floes.centroid
        pieces_list = StructArray{Floe{Float64}}(undef, 0)
        max_id = Subzero.timestep_ridging_rafting!(
            floes,
            pieces_list,
            3,
            periodic_domain,
            maximum(floes.id),
            ridge_settings,
            coupling_settings,
            simp_settings,
            consts,
            10,
        )
        # Make sure each piece is saved, ghosts are marked to remove and
        # don't interact with any other floes
        @test length(pieces_list) == 1
        @test total_mass == sum(floes.mass[1:3]) + sum(pieces_list.mass)
        # Make sure floe 1 parent ridged onto floe 3 
        @test h3 < floes.height[3]
        @test cent3 == floes.centroid[3]
        # Make sure floes 1 ghost ridged onto floe 2 and broke
        @test isapprox(h1, floes.height[1], atol = 1e-15) && h2 < floes.height[2]
        @test floes.height[1] == pieces_list.height[1]
        @test cent2 == floes.centroid[2]
        # Make sure IDs are correct
        @test max_id == 5
        @test floes.id == [4, 2, 3, 1] && pieces_list.id[1] == 5
        @test floes.parent_ids[1] == [1] && pieces_list.parent_ids[1] == [1]
        @test isempty(floes.parent_ids[2]) && isempty(floes.parent_ids[3])
        # Make sure ghost floe is marked for remove_floe_overlap
        @test floes.status[4].tag == Subzero.remove
    end
end