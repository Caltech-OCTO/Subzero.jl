@testset "Collisions" begin
    @testset "Normal Force" begin
        
    end
    @testset "Elastic Forces" begin
        
    end
    @testset "Frictional Forces" begin
        
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
        # flip border overlap from north to south and east to west and visa versa
        topo_arr = StructArray([TopographyElement(-coords1)])

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
        new_topo_arr = deepcopy(topo_arr)
        add_ghosts!(new_floe_arr, nonperiodic_domain)
        add_ghosts!(new_topo_arr, nonperiodic_domain)
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
        new_topo_arr = deepcopy(topo_arr)
        add_ghosts!(new_floe_arr, double_periodic_domain)
        add_ghosts!(new_topo_arr, double_periodic_domain)
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
        # Add ghosts for topography in both east-west and north-south directions
        @test -1e5 < new_topo_arr.centroid[1][1] < 1e5
        @test -1e5 < new_topo_arr.centroid[1][2] < 1e5
        @test new_topo_arr.coords[1] == Subzero.translate(topo_arr.coords[1], [2e5, 2e5])
        @test new_topo_arr.coords[2] == Subzero.translate(topo_arr.coords[1], [0, 2e5])
        @test new_topo_arr.coords[3] == Subzero.translate(topo_arr.coords[1], [2e5, 0])
        @test new_topo_arr.coords[4] == topo_arr.coords[1]
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