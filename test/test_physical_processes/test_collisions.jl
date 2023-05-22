@testset "Collisions" begin
    FT = Float64
    @testset "Floe-Floe Interactions" begin
        Δt = 10
        max_overlap = 0.55
        hmean = 0.25
        Δh = 0.0
        tri = Floe(
            FT,
            [[[0.0, 0.0], [1e4, 3e4], [2e4, 0], [0.0, 0.0]]],
            hmean,
            Δh,
        )
        tri.u = 0.1
        rect = Floe(
            FT,
            [[
                [0.0, 2.5e4],
                [0.0, 2.9e4],
                [2e4, 2.9e4],
                [2e4, 2.5e4],
                [0.0, 2.5e4],
            ]],
            hmean,
            Δh,
        )
        rect.v = -0.1
        cfloe = Floe(
            FT,
            [[
                [0.5e4, 2.7e4],
                [0.5e4, 3.5e4],
                [1.5e4, 3.5e4],
                [1.5e4, 2.7e4],
                [1.25e4, 2.7e4],
                [1.25e4, 3e4],
                [1e4, 3e4],
                [1e4, 2.7e4],
                [0.5e4, 2.7e4],
            ]],
            hmean,
            Δh,
        )
        cfloe.u = 0.3
        consts = Constants()
        # Triange tip intersected with a rectangle
        Subzero.floe_floe_interaction!(tri, 1, rect, 2, 2, consts, Δt, max_overlap)
        @test isapprox(tri.interactions[1, xforce], -64613382.47, atol = 1e-2)
        @test isapprox(tri.interactions[1, yforce], -521498991.51, atol = 1e-2)
        @test isapprox(tri.interactions[1, xpoint], 10000.00, atol = 1e-2)
        @test isapprox(tri.interactions[1, ypoint], 26555.55, atol = 1e-2)
        @test isapprox(tri.interactions[1, overlap], 8000000, atol = 1e-2)
        @test tri.status.tag != Subzero.fuse && rect.status.tag != Subzero.fuse
        @test isempty(tri.status.fuse_idx)
        Subzero.calc_torque!(tri)
        @test isapprox(tri.interactions[1, torque], 1069710443203.99, atol = 1e-2)
        tri.interactions = zeros(0, 7)
        # Sideways C intersected with rectangle so there are two areas of overlap
        Subzero.floe_floe_interaction!(cfloe, 1, rect, 2, 2, consts, Δt, max_overlap)
        @test isapprox(cfloe.interactions[1, xforce], -163013665.41, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, xforce], -81506832.70, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, yforce], 804819565.60, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, yforce], 402409782.80, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, xpoint], 7500.00, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, xpoint], 13750.00, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, ypoint], 28000.00, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, ypoint], 28000.00, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, overlap], 10000000, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, overlap], 5000000, atol = 1e-2)
        @test tri.status.tag != Subzero.fuse && rect.status.tag != Subzero.fuse
        Subzero.calc_torque!(cfloe)
        @test isapprox(cfloe.interactions[1, torque], -2439177121266.03, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, torque], 1295472581868.05, atol = 1e-2)
        cfloe.interactions = zeros(0, 7)
        # Floes overlapping more than 55%  - rectangle and shifted rectangle
        shift_rect = deepcopy(rect)
        shift_rect.coords = Subzero.translate(shift_rect.coords, 0.5e4, 0.0)
        Subzero.floe_floe_interaction!(rect, 1, shift_rect, 2, 2, consts, Δt, max_overlap)
        @test rect.status.tag == Subzero.fuse 
        @test rect.status.fuse_idx == [2]
        @test isempty(rect.interactions)
        empty!(rect.status.fuse_idx)
        rect.status.tag == Subzero.active
        # Floe 2 (small rectangle) is overlapping with rectangle by more than 55%
        small_rect = Floe(
            FT,
            [[
                [1.8e4, 2.7e4],
                [1.8e4, 2.8e4],
                [2.1e4, 2.8e4],
                [2.1e4, 2.7e4],
                [1.8e4, 2.7e4],
            ]],
            hmean,
            Δh,
        )
        Subzero.floe_floe_interaction!(
            rect,
            1,
            small_rect,
            2,
            2,
            consts,
            Δt,
            max_overlap,
        )
        @test rect.status.tag == Subzero.fuse 
        @test rect.status.fuse_idx == [2]
        
        # Overlapping barely floes - such a small overlap that forces are not calculated
        shift_rect = deepcopy(rect)
        shift_rect.coords = Subzero.translate(shift_rect.coords, 1.9999999e4, 0.0)
        Subzero.floe_floe_interaction!(shift_rect, 1, rect, 2, 2, consts, Δt, max_overlap)
        @test isempty(shift_rect.interactions)
    end


    @testset "Floe Boundary Interactions" begin
        Δt = 10
        max_overlap = 0.75
        Lx = 1e5
        Ly = Lx
        hmean = 0.25
        Δh = 0.0
        grid = RegRectilinearGrid(
            Float64,
            (-Lx, Lx),
            (-Ly, Ly),
            1e4,
            1e4,
        )
        nboundary = PeriodicBoundary(FT, North, grid)
        sboundary = PeriodicBoundary(FT, South, grid)
        eboundary = CollisionBoundary(FT, East, grid)
        wboundary = OpenBoundary(FT, West, grid)
        topo_elem = TopographyElement(
            FT,
            [[[1e4, 0.0], [0.0, 1e4], [1e4, 2e4], [2e4, 1e4], [1e4, 0.0]]],
        )
        domain = Domain(nboundary, sboundary, eboundary, wboundary, StructArray([topo_elem]))
        # Diagonal floe barely overlaping with eastern collision boundary
        efloe_small = Floe(
            FT,
            [[
                [9.5e4, 0.0],
                [9e4, 0.5e4],
                [10e4, 2.5e4],
                [10.05e4, 2e4],
                [9.5e4, 0.0],
            ]],
            hmean,
            Δh,
        )
        efloe_small.u = 0.5
        efloe_small.v = 0.25
        # Floe overlapping with eastern collision boundary by more than 75% to trigger overlap condition
        efloe_large = Floe(
            FT,
            [[
                [9e4, -7e4],
                [9e4, -5e4],
                [1.4e5, -5e4],
                [1.4e5, -7e4],
                [9e4, -7e4],
            ]],
            hmean,
            Δh,
        )
        efloe_large.u = 0.1
        efloe_large.v = -0.35
        # Floe overlapping with boundary with more than one region
        cfloe = Floe(
            FT,
            [[
                [9.5e4, 7e4],
                [9.5e4, 9e4],
                [1.05e5, 9e4],
                [1.05e5, 8.5e4],
                [9.9e4, 8.5e4],
                [9.9e4, 8e4],
                [1.05e5, 8e4],
                [1.05e5, 7e4],
                [9.5e4, 7e4],
            ]],
            hmean,
            Δh,
        )
        cfloe.v = -0.1
        # Floe overlapping with western open boundary
        wfloe = Floe(
            FT,
            [[
                [-9.75e4, 7e4],
                [-9.75e4, 5e4],
                [-10.05e4, 5e4],
                [-10.05e4, 7e4],
                [-9.75e4, 7e4],
            ]],
            hmean,
            Δh,
        )
        # Floe overlapping with northern periodic boundary
        nfloe = Floe(
            FT,
            [[
                [5e4, 9.75e4],
                [5e4, 10.05e4],
                [7e4, 10.05e4],
                [7e4, 9.75e4],
                [5e4, 9.75e4],
            ]],
            hmean,
            Δh,
        )
        # Floe overlapping with topography element
        tfloe = Floe(
            FT,
            [[
                [-0.5e4, 0.0],
                [-0.5e4, 0.75e4],
                [0.5e4, 0.75e4],
                [0.5e4, 0.0],
                [-0.5e4, 0.0],
            ]],
            hmean,
            Δh,
        )
        efloe_large.u = -0.4
        efloe_large.v = 0.2

        consts = Constants()

        # Test floe overlapping slightly with collision boundary
        Subzero.floe_domain_interaction!(efloe_small, domain, consts, Δt, max_overlap)
        @test efloe_small.interactions[1, floeidx] == Inf
        @test isapprox(efloe_small.interactions[1, xforce], -311304795.629, atol = 1e-3)
        @test isapprox(efloe_small.interactions[1, yforce], -23618874.648, atol = 1e-3)
        @test isapprox(efloe_small.interactions[1, overlap], 1704545.454, atol = 1e-3)
        @test isapprox(efloe_small.interactions[1, xpoint], 100166.666, atol = 1e-3)
        @test isapprox(efloe_small.interactions[1, ypoint], 21060.606, atol = 1e-3)

        Subzero.floe_domain_interaction!(cfloe, domain, consts, Δt, max_overlap)
        @test isapprox(cfloe.interactions[1, xforce], -2876118708.17, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, xforce], -5752237416.35, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, yforce], 575223741.63, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, yforce], 1150447483.27, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, xpoint], 102500, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, xpoint], 102500, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, ypoint], 87500, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, ypoint], 75000, atol = 1e-2)
        @test isapprox(cfloe.interactions[1, overlap], 25000000, atol = 1e-2)
        @test isapprox(cfloe.interactions[2, overlap], 50000000, atol = 1e-2)

        # Test floe overlapping >75% with collision boundary
        Subzero.floe_domain_interaction!(efloe_large, domain, consts, Δt, max_overlap)
        @test isempty(efloe_large.interactions)
        @test efloe_large.status.tag == Subzero.remove
        # Test floe overlapping >75% with collision boundary -> but max overlap is now 1
        Subzero.floe_domain_interaction!(efloe_large, domain, consts, Δt, 1.0)
        @test !isempty(efloe_large.interactions)
        # Test floe passing through open boundary is killed
        Subzero.floe_domain_interaction!(wfloe, domain, consts, Δt, max_overlap)
        @test wfloe.status.tag == Subzero.remove
        # Test floes not not interact with periodic boundary
        nfloe_copy = deepcopy(nfloe)
        Subzero.floe_domain_interaction!(nfloe, domain, consts, Δt, max_overlap)
        @test nfloe_copy.status.tag == nfloe.status.tag && nfloe_copy.interactions == nfloe.interactions
        # Test floe overlapping with topography -> different from Subzero since topography are now boundaries
        Subzero.floe_domain_interaction!(tfloe, domain, consts, Δt, max_overlap)
        @test tfloe.interactions[1, xforce] < 0
        @test tfloe.interactions[1, yforce] < 0

        # Test floe hitting more than one wall at once -> different from Subzero
        collision_domain = Domain(
            CollisionBoundary(FT, North, grid),
            CollisionBoundary(FT, South, grid),
            CollisionBoundary(FT, East, grid),
            CollisionBoundary(FT, West, grid),
        )
        corner_floe = Floe(
            FT,
            [[
                [9.5e4, 7e4],
                [9e4, 7.5e4],
                [10e4, 1.05e5],
                [10.05e4, 9.5e4],
                [9.5e4, 7e4],
            ]],
            hmean,
            Δh,
        )
        Subzero.floe_domain_interaction!(
            corner_floe,
            collision_domain,
            consts,
            Δt,
            max_overlap,
        )
        @test all(corner_floe.interactions[:, xforce] .<= 0)
        @test all(corner_floe.interactions[:, yforce] .<= 0)

    end
    
    @testset "Add Ghosts" begin
        FT = Float64
        Lx = 1e5
        grid = RegRectilinearGrid(
            FT,
            (-Lx, Lx),
            (-Lx, Lx),
            1e4,
            1e4,
        )
        nboundary = PeriodicBoundary(FT, North, grid)
        sboundary = PeriodicBoundary(FT, South, grid)
        eboundary = PeriodicBoundary(FT, East, grid)
        wboundary = PeriodicBoundary(FT, West, grid)

        # corner floe - overlaps with north and east boundary
        coords1 = [[[9.9e4, 9.9e4], [9.9e4, 1.02e5], [1.02e5, 1.02e5], [1.02e5, 9.9e4], [9.9e4, 9.9e4]]]
        # overlaps with western boundary
        coords2 = [[[-1.01e5, 7e4], [-1.01e5, 8e4], [-8e4, 8e4], [-8e4, 7e4], [-1.01e5, 7e4]]]
        # overlaps with northern boundary
        coords3 = [[[-2e4, 9.5e4], [-2e4, 1.1e5], [-1e4, 1.1e5], [-1e4, 9.5e4], [-2e4, 9.5e4]]]
        # doesn't overlap with any boundary
        coords4 = [[[0.0, 0.0], [0.0, 2e4], [2e4, 2e4], [2e4, 0.0], [0.0, 0.0]]]

        floe_arr = StructArray(
            [Floe(FT, c, 0.5, 0.0) for c in [coords1, coords2, coords3, coords4]],
        )
        for i in eachindex(floe_arr)
            floe_arr.id[i] = Float64(i)
        end

        nonperiodic_domain = Domain(
            OpenBoundary(FT, North, grid),
            OpenBoundary(FT, South, grid),
            OpenBoundary(FT, East, grid),
            OpenBoundary(FT, West, grid),
        )
        ew_periodic_domain = Domain(
            OpenBoundary(FT, North, grid),
            OpenBoundary(FT, South, grid),
            PeriodicBoundary(FT, East, grid),
            PeriodicBoundary(FT, West, grid),
        )
        ns_periodic_domain = Domain(
            PeriodicBoundary(FT, North, grid),
            PeriodicBoundary(FT, South, grid),
            OpenBoundary(FT, East, grid),
            OpenBoundary(FT, West, grid),
        )
        double_periodic_domain = Domain(
            PeriodicBoundary(FT, North, grid),
            PeriodicBoundary(FT, South, grid),
            PeriodicBoundary(FT, East, grid),
            PeriodicBoundary(FT, West, grid),
        )

        # Make sure nothing is added with non-periodic domain
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, nonperiodic_domain)
        @test new_floe_arr.coords == floe_arr.coords

        # Add ghost floes in east-west direction
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, ew_periodic_domain)
        @test -1e5 < new_floe_arr[1].centroid[1] < 1e5
        @test -1e5 < new_floe_arr[2].centroid[2] < 1e5
        @test new_floe_arr.coords[1] == Subzero.translate(floe_arr.coords[1], -2e5, 0.0)
        @test new_floe_arr.coords[2:4] == floe_arr.coords[2:4]
        @test new_floe_arr.coords[5] == floe_arr.coords[1]
        @test new_floe_arr.coords[6] == Subzero.translate(floe_arr.coords[2], 2e5, 0.0)
        @test new_floe_arr.id == [1, 2, 3, 4, 1, 2]
        @test new_floe_arr.ghost_id == [0, 0, 0, 0, 1, 1]
        @test new_floe_arr.ghosts[1] == [5]
        @test new_floe_arr.ghosts[2] == [6]
        @test new_floe_arr.ghosts[3:6] == [[], [], [], []]

        # Add ghost floes in the north-south direction
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, ns_periodic_domain)
        @test -1e5 < new_floe_arr[1].centroid[2] < 1e5
        @test -1e5 < new_floe_arr[3].centroid[2] < 1e5
        @test new_floe_arr.coords[1] == Subzero.translate(floe_arr.coords[1], 0.0, -2e5)
        @test new_floe_arr.coords[3] == Subzero.translate(floe_arr.coords[3], 0.0, -2e5)
        @test new_floe_arr.coords[[2, 4]] == floe_arr.coords[[2, 4]]
        @test new_floe_arr.coords[5] == floe_arr.coords[1]
        @test new_floe_arr.coords[6] == floe_arr.coords[3]
        @test new_floe_arr.id == [1, 2, 3, 4, 1, 3]
        @test new_floe_arr.ghost_id == [0, 0, 0, 0, 1, 1]
        @test new_floe_arr.ghosts[1] == [5]
        @test new_floe_arr.ghosts[3] == [6]
        @test new_floe_arr.ghosts[[2; 4:6]] == [[], [], [], []]

        # Add ghosts in both east-west and north-south directions
        new_floe_arr = deepcopy(floe_arr)
        add_ghosts!(new_floe_arr, double_periodic_domain)
        @test -1e5 < new_floe_arr.centroid[1][1] < 1e5
        @test -1e5 < new_floe_arr.centroid[1][2] < 1e5
        @test new_floe_arr.coords[1] == Subzero.translate(floe_arr.coords[1], -2e5, -2e5)
        @test new_floe_arr.coords[3] == Subzero.translate(floe_arr.coords[3], 0.0, -2e5)
        @test new_floe_arr.coords[[2, 4]] == floe_arr.coords[[2, 4]]
        @test new_floe_arr.coords[5] == floe_arr.coords[1]
        @test new_floe_arr.coords[6] == Subzero.translate(floe_arr.coords[2], 2e5, 0.0)
        @test new_floe_arr.coords[7] == Subzero.translate(floe_arr.coords[1], 0.0, -2e5)
        @test new_floe_arr.coords[8] == Subzero.translate(floe_arr.coords[1], -2e5, 0.0)
        @test new_floe_arr.coords[9] == floe_arr.coords[3]
        @test new_floe_arr.id == [1, 2, 3, 4, 1, 2, 1, 1, 3]
        @test new_floe_arr.ghost_id == [0, 0, 0, 0, 1, 1, 2, 3, 1]
        @test new_floe_arr.ghosts[1] == [5, 7, 8]
        @test new_floe_arr.ghosts[2] == [6]
        @test new_floe_arr.ghosts[3] == [9]
        @test new_floe_arr.ghosts[4:9] == [[], [], [], [], [], []]
    
    end
    @testset "Ghost Collisions" begin
        Δt = 10
        Lx = 1e5
        Ly = 1e5
        collision_settings = CollisionSettings()
        grid = RegRectilinearGrid(
            Float64,
            (-Lx, Lx),
            (-Lx, Lx),
            1e4,
            1e4,
        )
        double_periodic_domain = Domain(
            PeriodicBoundary(FT, North, grid),
            PeriodicBoundary(FT, South, grid),
            PeriodicBoundary(FT, East, grid),
            PeriodicBoundary(FT, West, grid),
        )
        # Parent-parent collison (parents are touching)
        coords1 = splitdims([Lx/2 Lx/2 3*Lx/4 3*Lx/4 Lx+10000 Lx+10000; Ly/2 Ly+10000 Ly+10000 3*Ly/4 3*Ly/4 Ly/2])
        th = 0:pi/50:2*pi
        r = Ly/4+1000
        coords2 = invert([r * cos.(th) .+ (Lx-1), r * sin.(th) .+ (Ly-1)])

        floe_arr = StructArray(Floe(FT, [c], 0.5, 0.0) for c in [coords1, coords2])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = Float64(i)
        end
        spinlock = Threads.SpinLock()
        Subzero.timestep_collisions!(
            floe_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        xforce_vals = abs(floe_arr[1].collision_force[1])
        yforce_vals = abs(floe_arr[1].collision_force[2])
        f1_torque = floe_arr[1].collision_trq
        f2_torque = floe_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(
            floe_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        # 1 and 2 are the "parent" floes - floe 1 and floe 2 interact
        @test xforce_vals == abs(floe_arr[1].collision_force[1]) == abs(floe_arr[2].collision_force[1])
        @test yforce_vals == abs(floe_arr[2].collision_force[2]) == abs(floe_arr[2].collision_force[2])
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq
        # All other ghost floes aren't calculated
        @test first.(floe_arr[3:end].collision_force) == zeros(6)
        @test last.(floe_arr[3:end].collision_force) == zeros(6)
        @test floe_arr[3:end].collision_trq == zeros(6)

        # Ghost-Ghost collision (parents aren't touching, only ghosts touch)
        coords1 = splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))
        coords2 = splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], -[7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))
        floe_arr = StructArray(Floe(FT, [c], 0.5, 0.0) for c in [coords1, coords2])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = i
        end
        trans_arr = StructArray([
            Floe(
                FT,
                Subzero.translate([coords1],
                0.0, -2Ly),
                0.5,
                0.0,
            ),
            Floe(
                FT,
                Subzero.translate([coords2], 2Lx, 0.0),
                0.5,
                0.0,
            ),
        ])
        for i in eachindex(trans_arr)
            trans_arr.id[i] = i
        end
        Subzero.timestep_collisions!(
            trans_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        xforce_vals = abs(trans_arr[1].collision_force[1])
        yforce_vals = abs(trans_arr[1].collision_force[2])
        f1_torque = trans_arr[1].collision_trq
        f2_torque = trans_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(
            floe_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        # floes 1 and 2 are the parents - floe 4 is floe 1's ghost and floe 3 is
        # floe 2's ghost - floe 3 and 4 collide 
        @test repeat([xforce_vals], 2) == abs.(first.(floe_arr.collision_force[1:2]))
        @test repeat([yforce_vals], 2) == abs.(last.(floe_arr.collision_force[1:2]))
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq
        # interactions copied from ghosts
        @test floe_arr[1].interactions[:, [1:5; 7]] == floe_arr[4].interactions[:, [1:5; 7]]
        @test floe_arr[2].interactions[:, [1:5; 7]] == floe_arr[3].interactions[:, [1:5; 7]]

        # Parent-Ghost Collision
        coords1 = splitdims(vcat([5*Lx/8 5*Lx/8 3*Lx/4 3*Lx/4].+1000, [3*Ly/4 5*Ly/4 5*Ly/4 3*Ly/4]))
        coords2 = splitdims(vcat(-[5*Lx/4 5*Lx/4 3*Lx/4-1000 3*Lx/4-1000], [7*Lx/8 3*Lx/4-1000 3*Lx/4-1000 7*Lx/8]))
        floe_arr = StructArray(Floe(FT, [c], 0.5, 0.0) for c in [coords1, coords2])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = Float64(i)
        end
        trans_arr = StructArray([Floe(FT, Subzero.translate([coords1], -2Lx, 0.0), 0.5, 0.0),
                                 Floe(FT, [coords2], 0.5, 0.0)])
        for i in eachindex(trans_arr)
            trans_arr.id[i] = Float64(i)
        end
        Subzero.timestep_collisions!(
            trans_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        xforce_vals = abs(trans_arr[1].collision_force[1])
        yforce_vals = abs(trans_arr[1].collision_force[2])
        f1_torque = trans_arr[1].collision_trq
        f2_torque = trans_arr[2].collision_trq
        add_ghosts!(floe_arr, double_periodic_domain)
        Subzero.timestep_collisions!(
            floe_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        # Floe 1's ghost if floe 4 and floe 2's ghost is floe 3 and floe 3 and floe 1 interact
        @test repeat([xforce_vals], 2) == abs.(first.(floe_arr.collision_force[1:2]))
        @test repeat([yforce_vals], 2) == abs.(last.(floe_arr.collision_force[1:2]))
        @test f1_torque == floe_arr[1].collision_trq
        @test f2_torque == floe_arr[2].collision_trq
        @test floe_arr[2].interactions[:, [1:5; 7]] == floe_arr[3].interactions[:, [1:5; 7]]
        @test isempty(floe_arr[4].interactions)

        # Parent and ghosts hitting the same floe
        # small rectangle in the corner that has 3 ghosts in all other corners
        small_rect = Floe(
            FT,
            [[
                [-1.1e5, -1.1e5],
                [-1.1e5, -9.5e4],
                [-9.5e4, -9.5e4],
                [-9.5e4, -1.1e5],
                [-1.1e5, -1.1e5],
            ]],
            0.5,
            0.0,
        )
        # triangle in the middle of the domain with no ghosts - touches 3/4 corners
        large_tri = Floe(
            FT,
            [[
                [-1e5, -1e5],
                [-1e5, 1e5],
                [1e5, -1e5],
                [-1e5, -1e5],
            ]],
            0.5,
            0.0,
        )
        # rectangle along south boundary, ghost along north boundary
        bound_rect =  Floe(
            FT, 
            [[
                [-9.8e4, -1.1e5],
                [-9.8e4, -9.5e4],
                [9.8e4, -9.5e4],
                [9.8e4, -1.1e5],
                [-9.8e4, -1.1e5],
            ]],
            0.5,
            0.0,
        )
        floe_arr = StructArray([small_rect, large_tri])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = Float64(i)
        end
        add_ghosts!(floe_arr, double_periodic_domain)
        @test length(floe_arr) == 5
        Subzero.timestep_collisions!(
            floe_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        @test size(floe_arr[1].interactions)[1] == 3
        @test size(floe_arr[2].interactions)[1] == 3
        @test floe_arr[1].interactions[1, Subzero.xforce] != floe_arr[1].interactions[2, Subzero.xforce] && floe_arr[1].interactions[1, Subzero.xforce] != floe_arr[1].interactions[3,Subzero.xforce]
        @test floe_arr[1].interactions[1, Subzero.yforce] != floe_arr[1].interactions[2, Subzero.yforce] && floe_arr[1].interactions[1, Subzero.yforce] != floe_arr[1].interactions[3, Subzero.yforce]

        floe_arr = StructArray([small_rect, bound_rect])
        for i in eachindex(floe_arr)
            floe_arr.id[i] = Float64(i)
        end
        add_ghosts!(floe_arr, double_periodic_domain)
        @test length(floe_arr) == 6
        Subzero.timestep_collisions!(
            floe_arr,
            2,
            double_periodic_domain,
            Subzero.Constants(),
            Δt,
            collision_settings,
            spinlock,
        )
        @test size(floe_arr[1].interactions)[1] == 2
        @test size(floe_arr[2].interactions)[1] == 2
        @test floe_arr[1].interactions[1, Subzero.xpoint] != floe_arr[1].interactions[2, Subzero.xpoint]
        @test floe_arr[1].interactions[1, Subzero.ypoint] == floe_arr[1].interactions[2, Subzero.ypoint]
    end
end