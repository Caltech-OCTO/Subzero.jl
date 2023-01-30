@testset "Collisions" begin
    @testset "Floe-Floe Interactions" begin
        
    end


    @testset "Floe Boundary Interactions" begin
        Lx = 1e5
        Ly = Lx
        h_mean = 0.25
        Δh = 0.0
        grid = RegRectilinearGrid(-Lx, Lx, -Ly, Ly, 1e4, 1e4)
        nboundary = PeriodicBoundary(grid, North())
        sboundary = PeriodicBoundary(grid, South())
        eboundary = CollisionBoundary(grid, East())
        wboundary = OpenBoundary(grid, West())
        topo_elem = TopographyElement([[[1e4, 0.0], [0.0, 1e4], [1e4, 2e4], [2e4, 1e4], [1e4, 0.0]]])
        domain = Domain(nboundary, sboundary, eboundary, wboundary, StructArray([topo_elem]))
        # Diagonal floe barely overlaping with eastern collision boundary
        efloe_small = Floe([[[9.5e4, 0.0], [9e4, 0.5e4], [10e4, 2.5e4], [10.05e4, 2e4], [9.5e4, 0.0]]], h_mean, Δh)
        efloe_small.u = 0.5
        efloe_small.v = 0.25
        # Floe overlapping with eastern collision boundary by more than 75% to trigger overlap condition
        efloe_large = Floe([[[9e4, -7e4], [9e4, -5e4], [1.4e5, -5e4], [1.4e5, -7e4], [9e4, -7e4]]], h_mean, Δh)
        efloe_large.u = 0.1
        efloe_large.v = -0.35
        # Floe overlapping with western open boundary
        wfloe = Floe([[[-9.75e4, 7e4], [-9.75e4, 5e4], [-10.05e4, 5e4], [-10.05e4, 7e4], [-9.75e4, 7e4]]], h_mean, Δh)
        # Floe overlapping with northern periodic boundary
        nfloe = Floe([[[5e4, 9.75e4], [5e4, 10.05e4], [7e4, 10.05e4], [7e4, 9.75e4], [5e4, 9.75e4]]], h_mean, Δh)
        # Floe overlapping with topography element
        tfloe = Floe([[[-0.5e4, 0.0], [-0.5e4, 0.75e4], [0.5e4, 0.75e4], [0.5e4, 0.0], [-0.5e4, 0.0]]], h_mean, Δh)
        efloe_large.u = -0.4
        efloe_large.v = 0.2
        floe_arr = StructArray([efloe_small, efloe_large, wfloe, nfloe, tfloe])
        consts = Constants()

        # Test floe overlapping slightly with collision boundary
        Subzero.floe_domain_interaction!(efloe_small, domain, consts, 10)
        @test efloe_small.interactions[2, "floeidx"] == Inf
        @test isapprox(efloe_small.interactions[2, "xforce"], -311304795.629, atol = 1e-3)
        # if we want to zero out then should be 0
        @test isapprox(efloe_small.interactions[2, "yforce"], -23618874.648, atol = 1e-3)
        @test isapprox(efloe_small.interactions[2, "overlap"], 1704545.454, atol = 1e-3)
        @test isapprox(efloe_small.interactions[2, "xpoint"], 100166.666, atol = 1e-3)
        @test isapprox(efloe_small.interactions[2, "ypoint"], 21060.606, atol = 1e-3)

        # Test floe overlapping >75% with collision boundary
        Subzero.floe_domain_interaction!(efloe_large, domain, consts, 10)
        @test size(efloe_large.interactions, 1) == 1  # No interactions
        @test efloe_large.alive == 0
        # Test floe passing through open boundary is killed
        Subzero.floe_domain_interaction!(wfloe, domain, consts, 10)
        @test wfloe.alive == 0
        # Test floes not not interact with periodic boundary
        nfloe_copy = deepcopy(nfloe)
        Subzero.floe_domain_interaction!(nfloe, domain, consts, 10)
        @test nfloe_copy.alive == nfloe.alive &&  nfloe_copy.interactions == nfloe.interactions
        # Test floe overlapping with topography
        Subzero.floe_domain_interaction!(tfloe, domain, consts, 10)
        # TODO: this has issues. Need to check on it. 

    end
    
    @testset "Add Ghosts" begin
        Lx = 1e5
        grid = RegRectilinearGrid(-Lx, Lx, -Lx, Lx, 1e4, 1e4)
        nboundary = PeriodicBoundary(grid, North())
        sboundary = PeriodicBoundary(grid, South())
        eboundary = PeriodicBoundary(grid, East())
        wboundary = PeriodicBoundary(grid, West())

        # corner floe - overlaps with north and east boundary
        coords1 = [[[9.9e4, 9.9e4], [9.9e4, 1.02e5], [1.02e5, 1.02e5], [1.02e5, 9.9e4], [9.9e4, 9.9e4]]]
        # overlaps with western boundary
        coords2 = [[[-1.01e5, 7e4], [-1.01e5, 8e4], [-8e4, 8e4], [-8e4, 7e4], [-1.01e5, 7e4]]]
        # overlaps with northern boundary
        coords3 = [[[-2e4, 9.5e4], [-2e4, 1.1e5], [-1e4, 1.1e5], [-1e4, 9.5e4], [-2e4, 9.5e4]]]
        # doesn't overlap with any boundary
        coords4 = [[[0.0, 0.0], [0.0, 2e4], [2e4, 2e4], [2e4, 0.0], [0.0, 0.0]]]

        floe_arr = StructArray([Floe(c, 0.5, 0.0) for c in [coords1, coords2, coords3, coords4]])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = i
        end

        nonperiodic_domain = Domain(OpenBoundary(grid, North()), OpenBoundary(grid, South()),
                                    OpenBoundary(grid, East()), OpenBoundary(grid, West()), topo_arr)

        ew_periodic_domain = Domain(OpenBoundary(grid, North()), OpenBoundary(grid, South()),
                                    PeriodicBoundary(grid, East()), PeriodicBoundary(grid, West()), topo_arr)

        ns_periodic_domain = Domain(PeriodicBoundary(grid, North()), PeriodicBoundary(grid, South()),
                                    OpenBoundary(grid, East()), OpenBoundary(grid, West()), topo_arr)

        double_periodic_domain = Domain(PeriodicBoundary(grid, North()), PeriodicBoundary(grid, South()),
                                        PeriodicBoundary(grid, East()), PeriodicBoundary(grid, West()), topo_arr)

        # Make sure nothing is added with non-periodic domain
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, nonperiodic_domain)
        @test new_floe_arr.coords == floe_arr.coords
        @test new_topo_arr.coords == topo_arr.coords

        # Add ghost floes in east-west direction
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, ew_periodic_domain)
        @test -1e5 < new_floe_arr[1].centroid[1] < 1e5
        @test -1e5 < new_floe_arr[2].centroid[2] < 1e5
        @test new_floe_arr.coords[1] == Subzero.translate(floe_arr.coords[1], [-2e5, 0.0])
        @test new_floe_arr.coords[2:4] == floe_arr.coords[2:4]
        @test new_floe_arr.coords[5] == floe_arr.coords[1]
        @test new_floe_arr.coords[6] == Subzero.translate(floe_arr.coords[2], [2e5, 0.0])
        @test new_floe_arr.id == [1, 2, 3, 4, -1, -2]
        @test new_floe_arr.ghosts[1] == [5]
        @test new_floe_arr.ghosts[2] == [6]
        @test new_floe_arr.ghosts[3:6] == [[], [], [], []]

        # Add ghost floes in the north-south direction
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, ns_periodic_domain)
        @test -1e5 < new_floe_arr[1].centroid[2] < 1e5
        @test -1e5 < new_floe_arr[3].centroid[2] < 1e5
        @test new_floe_arr.coords[1] == Subzero.translate(floe_arr.coords[1], [0.0, -2e5])
        @test new_floe_arr.coords[3] == Subzero.translate(floe_arr.coords[3], [0.0, -2e5])
        @test new_floe_arr.coords[[2, 4]] == floe_arr.coords[[2, 4]]
        @test new_floe_arr.coords[5] == floe_arr.coords[1]
        @test new_floe_arr.coords[6] == floe_arr.coords[3]
        @test new_floe_arr.id == [1, 2, 3, 4, -1, -3]
        @test new_floe_arr.ghosts[1] == [5]
        @test new_floe_arr.ghosts[3] == [6]
        @test new_floe_arr.ghosts[[2; 4:6]] == [[], [], [], []]

        # Add ghosts in both east-west and north-south directions
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, double_periodic_domain)
        @test -1e5 < new_floe_arr.centroid[1][1] < 1e5
        @test -1e5 < new_floe_arr.centroid[1][2] < 1e5
        @test new_floe_arr.coords[1] == Subzero.translate(floe_arr.coords[1], [-2e5, -2e5])
        @test new_floe_arr.coords[3] == Subzero.translate(floe_arr.coords[3], [0.0, -2e5])
        @test new_floe_arr.coords[[2, 4]] == floe_arr.coords[[2, 4]]
        @test new_floe_arr.coords[5] == floe_arr.coords[1]
        @test new_floe_arr.coords[6] == Subzero.translate(floe_arr.coords[2], [2e5, 0.0])
        @test new_floe_arr.coords[7] == Subzero.translate(floe_arr.coords[1], [0.0, -2e5])
        @test new_floe_arr.coords[8] == Subzero.translate(floe_arr.coords[1], [-2e5, 0.0])
        @test new_floe_arr.coords[9] == floe_arr.coords[3]
        @test new_floe_arr.id == [1, 2, 3, 4, -1, -2, -1, -1, -3]
        @test new_floe_arr.ghosts[1] == [5, 7, 8]
        @test new_floe_arr.ghosts[2] == [6]
        @test new_floe_arr.ghosts[3] == [9]
        @test new_floe_arr.ghosts[4:9] == [[], [], [], [], [], []]
    
    end
    @testset "Ghost Collisions" begin
        Lx = 1e5
        Ly = 1e5
        grid = RegRectilinearGrid(-Lx, Lx, -Lx, Lx, 1e4, 1e4)
        double_periodic_domain = Domain(PeriodicBoundary(grid, North()), PeriodicBoundary(grid, South()),
                                        PeriodicBoundary(grid, East()), PeriodicBoundary(grid, West()))
        # Parent-parent collison (parents are touching)
        coords1 = splitdims([Lx/2 Lx/2 3*Lx/4 3*Lx/4 Lx+10000 Lx+10000; Ly/2 Ly+10000 Ly+10000 3*Ly/4 3*Ly/4 Ly/2])
        th = 0:pi/50:2*pi
        r = Ly/4+1000
        coords2 = invert([r * cos.(th) .+ (Lx-1), r * sin.(th) .+ (Ly-1)])

        floe_arr = StructArray(Floe([c], 0.5, 0.0) for c in [coords1, coords2])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = i
        end
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
        xforce = abs(floe_arr[1].collision_force[1])
        yforce = abs(floe_arr[1].collision_force[2])
        f1_torque = floe_arr[1].collision_trq
        f2_torque = floe_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
        # 1 and 2 are the "parent" floes - floe1 interacts with floe 2's ghost floe (floe 4)
        @test xforce == abs(floe_arr[1].collision_force[1]) == abs(floe_arr[2].collision_force[1])
        @test yforce == abs(floe_arr[2].collision_force[2]) == abs(floe_arr[2].collision_force[2])
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq
        # All other ghost floes aren't calculated
        @test first.(floe_arr[3:end].collision_force) == zeros(6)
        @test last.(floe_arr[3:end].collision_force) == zeros(6)
        @test floe_arr[3:end].collision_trq == zeros(6)

        # Ghost-Ghost collision (parents aren't touching, only ghosts touch)
        coords1 = splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))
        coords2 = splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], -[7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))
        floe_arr = StructArray(Floe([c], 0.5, 0.0) for c in [coords1, coords2])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = i
        end
        trans_arr = StructArray([Floe(Subzero.translate([coords1], [0.0, -2Ly]), 0.5, 0.0),
                                 Floe(Subzero.translate([coords2], [2Lx, 0.0]), 0.5, 0.0)])
        for i in eachindex(trans_arr)
            trans_arr.id[i] = i
        end
        Subzero.timestep_collisions!(trans_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
        xforce = abs(trans_arr[1].collision_force[1])
        yforce = abs(trans_arr[1].collision_force[2])
        f1_torque = trans_arr[1].collision_trq
        f2_torque = trans_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
        @test repeat([xforce], 4) == abs.(first.(floe_arr.collision_force[1:4]))
        @test repeat([yforce], 4) == abs.(last.(floe_arr.collision_force[1:4]))
        @test f1_torque == floe_arr[1].collision_trq == floe_arr[4].collision_trq
        @test f2_torque == floe_arr[2].collision_trq == floe_arr[3].collision_trq
        @test size(floe_arr[1].interactions, 1) == 1
        @test size(floe_arr[2].interactions, 1) == 1

        # Parent-Ghost Collision
        coords1 = splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))
        coords2 = splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], [7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))
        floe_arr = StructArray(Floe([c], 0.5, 0.0) for c in [coords1, coords2])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = i
        end
        trans_arr = StructArray([Floe(Subzero.translate([coords1], [-2Lx, 0.0]), 0.5, 0.0),
                                 Floe([coords2], 0.5, 0.0)])
        for i in eachindex(trans_arr)
            trans_arr.id[i] = i
        end
        Subzero.timestep_collisions!(trans_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
        xforce = abs(trans_arr[1].collision_force[1])
        yforce = abs(trans_arr[1].collision_force[2])
        f1_torque = trans_arr[1].collision_trq
        f2_torque = trans_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(floe_arr, 2, double_periodic_domain, zeros(Int, 2), zeros(Int, 2), Subzero.Constants(), 10)
        @test repeat([xforce], 3) == abs.(first.(floe_arr.collision_force[1:3]))
        @test repeat([yforce], 3) == abs.(last.(floe_arr.collision_force[1:3]))
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq == floe_arr[3].collision_trq
        @test size(floe_arr[2].interactions, 1) == 1
        @test size(floe_arr[4].interactions, 1) == 1
    end

end