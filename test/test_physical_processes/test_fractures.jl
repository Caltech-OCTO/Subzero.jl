@testset "Fractures" begin
    @testset "Fracture Criteria" begin
        # Test NoFracturee criteria
        @test NoFracture() isa NoFracture
        # Test HiblerYieldCurve criteria
        @test_throws ArgumentError HiblerYieldCurve(2.25, 20.0, [[[0.0, 0.0]]])
        @test HiblerYieldCurve(
            2.25,
            20.0,
            [[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]]
        ) isa HiblerYieldCurve
        # Test calculate_hibler
        hibler_verts = Subzero.calculate_hibler(0.5, 5e5, -1)
        hibler_poly = LG.Polygon(hibler_verts)
        @test isapprox(LG.area(hibler_poly), 49054437859.374, atol = -1e3)
        @test all(isapprox.(
            Subzero.find_poly_centroid(hibler_poly),
            [-1.25e5, -1.25e5],
            atol = 1e-3
        ))
        x_verts, y_verts = Subzero.separate_xy(hibler_verts)
        @test all(isapprox.(
            extrema(x_verts),
            [-264743.588, 14727.999],
            atol = 1e-3
        ))
        @test all(isapprox.(
            extrema(y_verts),
            [-264743.588, 14727.999],
            atol = 1e-3
        ))
        hibler_verts = Subzero.calculate_hibler(0.25, 2.25e5, 20.0)
        hibler_poly = LG.Polygon(hibler_verts)
        @test isapprox(LG.area(hibler_poly), 2483380916.630, atol = -1e3)
        @test all(isapprox.(
            Subzero.find_poly_centroid(hibler_poly),
            [-28125, -28125],
            atol = 1e-3
        ))
        x_verts, y_verts = Subzero.separate_xy(hibler_verts)
        @test all(isapprox.(
            extrema(x_verts),
            [-59567.307, 3313.799],
            atol = 1e-3
        ))
        @test all(isapprox.(
            extrema(y_verts),
            [-59567.307, 3313.799],
            atol = 1e-3
        ))
        # Test Mohr's Cone
        @test typeof(Subzero.MohrsCone(Float32)) <: MohrsCone{Float32}
        @test typeof(Subzero.MohrsCone(Float64)) <: MohrsCone{Float64}

        # Float64 Mohr's Cone with q, σc, σ11
        mohrs_verts_64 = Subzero.calculate_mohrs(5.2, 2.5e5, -3.375e4)
        @test all(isapprox.(
            mohrs_verts_64[1],
            [
                [59523.809, 59523.809],
                [33750.0, -74500.0],
                [-74500.0, 33750.0],
                [59523.809, 59523.809],
            ],
            atol = 1e-3
        ))
        # Float32 Mohr's Cone with q, σc, σ11
        mohrs_verts_32 = Subzero.calculate_mohrs(5.2, 2.5e5, 1.5e5)
        @test all(isapprox.(
            mohrs_verts_32[1],
            [
                [59523.809, 59523.809],
                [-150000.0, -1.03e6],
                [-1.03e6, -150000.0],
                [59523.809, 59523.809],
            ],
            atol = 1e-3,
        ))
        # Float64 Mohr's Cone with σ1, σ2, σ11, σ22
        mohrs_verts_coords = Subzero.calculate_mohrs(
            5.95e4,
            5.95e4,
            -1.5e5,
            -1.03e6,
        )
        # Test update criteria for Hibler's curve
        floes = StructArray([Floe(
            [[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]],
            0.25,  # Floe has a height of 0.25
            0.0,
        )])
        yield_curve = HiblerYieldCurve(floes)
        verts = deepcopy(yield_curve.vertices)
        @test yield_curve isa HiblerYieldCurve
        @test yield_curve.pstar == 2.25e5 && yield_curve.c == 20
        floes.height .= 0.5
        Subzero.update_criteria!(yield_curve, floes)
        @test verts != yield_curve.vertices
        # Test update criteria for Mohr's cone
        cone_curve = MohrsCone()
        verts = cone_curve.vertices
        Subzero.update_criteria!(cone_curve, floes)
        @test verts == cone_curve.vertices
    end
    @testset "Fractures Floes" begin
        # Fracture tests depend on these floes and settings
        frac_stress = [-29955.396 -3428.008; -3428.008	-1942.0464]
        frac_deform_floe = Floe(
            [[
                [-50548.186, -49995.968],
                [-50550.745, -37790.078],
                [-20856.010, -32518.566],
                [-20929.577, -49989.757],
                [-50548.186, -49995.968],
            ]],
            0.25,
            0.0,
            u = 0.1,
            v = -0.2,
            ξ = 0.05,
        )
        frac_floe = deepcopy(frac_deform_floe)  # Without interactions, won't deform
        no_frac_floe = Floe(  # This floe is colliding with frac_deform_floe
            [[
                [1467.795, -25319.563],
                [1664.270, -25640.216],
                [-1105.179, -33458.936],
                [-17529.019, -50035.583],
                [-21193.828, -50088.777],
                [-21370.170, -32618.322],
                [-21247.656, -31077.536],
                [-12818.593, -27031.048],
                [1467.795, -25319.563],
            ]],
            0.25,
            0.0,
        )
        no_frac_small = Floe(  # This floe is too small to fracture or deform
            [[
                [1e3, 1e3],
                [1e3, 1.5e3],
                [1.5e3, 1.5e3],
                [1.5e3, 1e3],
                [1e3, 1e3],
            ]],
            0.25,
            0.0,
        )
        frac_deform_floe.stress = frac_stress
        frac_deform_floe.interactions = collect([
            3,
            -279441968.984,
            -54223517.438,
            -21091.0918258529,
            -40358.0042297616,
            -148920620521.112,
            6795329.38154967,
        ]')
        frac_deform_floe.num_inters = 1
        frac_deform_floe.p_dudt = 0.11
        frac_floe.stress = frac_stress
        no_frac_small.stress = frac_stress

        floes = StructArray([
            frac_deform_floe, frac_floe, no_frac_floe, no_frac_small
        ])
        floes.id .= collect(1:4)
        frac_settings = FractureSettings(
            fractures_on = true,
            criteria = HiblerYieldCurve(floes),
            Δt = 75,
            deform_on = true,
        )

        # Test determine_fractures
        frac_idx = Subzero.determine_fractures(
             floes,
             HiblerYieldCurve(floes),
             1e6
        )
        # First floe fractures, second is too small, third stress is too small
        @test frac_idx == [1, 2]
        
        # Test deform_floe!
        floe1_copy = deepcopy(floes[1])
        colliding_coords = no_frac_floe.coords
        deforming_forces = frac_deform_floe.interactions[xforce:yforce]
        init_overlap = LG.area(LG.intersection(
            LG.Polygon(floe1_copy.coords),
            LG.Polygon(colliding_coords),
        ))
        Subzero.deform_floe!(
            floe1_copy,
            colliding_coords,
            deforming_forces,
            CouplingSettings(),
            Constants(),
            10,
            Xoshiro(1),
        )
        @test init_overlap > LG.area(LG.intersection(
            LG.Polygon(floe1_copy.coords),  # These coords have changed
            LG.Polygon(colliding_coords),
        ))
        
        @test all(isapprox.( 
            floe1_copy.centroid,
            [-35115.567, -42531.500],
            atol = 2e-1))
        @test isapprox(floe1_copy.area, 431454521, atol = 10)

        # Test split_floe
        new_floes = Subzero.split_floe(
            floes[1],
            Xoshiro(3),
            frac_settings,
            CouplingSettings(),
            Constants(),
            10,
        ) 
        # Test that the pieces all fit within original floe
        og_floe_poly = LG.Polygon(floes.coords[1])
        new_floes_polys = LG.MultiPolygon(new_floes.coords)
        @test isapprox(
            LG.area(LG.intersection(new_floes_polys, og_floe_poly)),
            LG.area(og_floe_poly),
            atol = 1e-6,
        )
        # Conserve mass
        @test isapprox(sum(new_floes.mass), floes.mass[1], atol = 1e-4)
        # Linear velocities unchanged -> linear momentum conserved
        @test all(new_floes.u .== floes.u[1])
        @test all(new_floes.v .== floes.v[1])
        @test all(new_floes.p_dxdt .== floes.p_dxdt[1])
        @test all(new_floes.p_dudt .== floes.p_dudt[1])
        @test all([all(f.strain .== floes.strain[1]) for f in new_floes])
        # Angular momentum is conserved
        @test isapprox(
            sum(new_floes.ξ .* new_floes.moment),
            floes.ξ[1] * floes.moment[1],
            atol = 1e-8,
        )
        @test isapprox(
            sum(new_floes.p_dαdt .* new_floes.moment),
            floes.p_dαdt[1] * floes.moment[1],
            atol = 1e-8,
        )

        # Test fracture_floes!
        max_idx = Subzero.fracture_floes!(
            floes,
            4,  # start with 4 floes
            Xoshiro(3),
            frac_settings,
            CouplingSettings(),
            SimplificationSettings(),
            Constants(),
            10,
        )
        @test max_idx == 10
        @test length(floes) == 8
        @test all(floes.id .== [3, 4, 5, 6, 7, 8, 9, 10])
        @test all(floes.parent_ids .== [[], [], [2], [2], [2], [1], [1], [1]])
    end
end