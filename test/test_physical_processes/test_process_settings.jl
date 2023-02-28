@testset "ProcessInfo" begin
    @testset "CouplingSettings" begin
        # Test default
        default_info = CouplingSettings()
        @test default_info.coupling_on &&
            default_info.Δt == 10 &&
            !default_info.calc_ocnτ_on
        # Test custom correct
        custom_info = CouplingSettings(
            coupling_on = false,
            Δt = 20,
            calc_ocnτ_on = false,
        )
        @test !custom_info.coupling_on &&
            custom_info.Δt == 20 &&
            !custom_info.calc_ocnτ_on
        # Test negative timestep -> turns coupling off so turn calc_ocnτ off 
        warn_str1 = "Coupling can't occur on a multiple of negative timesteps. Turning coupling off."
        warn_str2 = "Can't calculate stresses on ocean from ice and atmosphere without coupling. Turning calc_ocnτ_on off."
        neg_Δt_info = @test_logs (:warn, warn_str1) (:warn, warn_str2) CouplingSettings(
                coupling_on = true,
                Δt = -20,
                calc_ocnτ_on = true,
            )
        @test !neg_Δt_info.coupling_on &&
            neg_Δt_info.Δt == -20 &&
            !neg_Δt_info.calc_ocnτ_on
    end

    @testset "CollisionSettings" begin
        # Test default
        default_info = CollisionSettings()
        @test default_info.collisions_on &&
            default_info.floe_floe_max_overlap == 0.55 &&
            default_info.floe_domain_max_overlap == 0.75
        # Test custom false
        custom_info = CollisionSettings(
            collisions_on = false,
            floe_floe_max_overlap = -0.1,
            floe_domain_max_overlap = 2.0
        )
        @test !custom_info.collisions_on &&
            custom_info.floe_floe_max_overlap == -0.1 &&
            custom_info.floe_domain_max_overlap == 2
        # Test overlaps out of bounds
        f2large_str = "The maximum collisin overlap between floes can't be greater than 1. Setting to 1."
        f2small_str = "The maximum collisin overlap between floes can't be less than 0. Setting to 0."
        d2large_str = "The maximum collisin overlap between floes and the domain can't be greater than 1. Setting to 1."
        d2small_str = "The maximum collisin overlap between floes and the domain can't be less than 0. Setting to 0."
        out_of_bounds = @test_logs (:warn, f2large_str) (:warn, d2small_str) CollisionSettings(
                collisions_on = true,
                floe_floe_max_overlap = 1.2,
                floe_domain_max_overlap = -0.1,
            )
        @test out_of_bounds.collisions_on &&
            out_of_bounds.floe_floe_max_overlap == 1.0 &&
            out_of_bounds.floe_domain_max_overlap == 0.0
        
        out_of_bounds = @test_logs (:warn, f2small_str) (:warn, d2large_str) CollisionSettings(
                collisions_on = true,
                floe_floe_max_overlap = -0.1,
                floe_domain_max_overlap = 1.2,
            )
        @test out_of_bounds.collisions_on &&
            out_of_bounds.floe_floe_max_overlap == 0.0 &&
            out_of_bounds.floe_domain_max_overlap == 1.0

    end
    @testset "FractureSettings" begin
        # Test default
        default_info = FractureSettings()
        @test !default_info.fractures_on &&
            default_info.criteria isa NoFracture &&
            default_info.Δt == 0 &&
            !default_info.deform_on &&
            default_info.npieces == 3
        # Test custom correct
        criteria = HiblerYieldCurve(
            0.0,
            0.0,
            [[[0.0, 0.0], [0, 1], [1 ,1], [1, 0]]],
        )
        custom_info = FractureSettings(
            fractures_on = true,
            criteria = criteria,
            Δt = 200,
            deform_on = true,
            npieces = 4
        )
        @test custom_info.fractures_on &&
            custom_info.criteria isa HiblerYieldCurve &&
            custom_info.Δt == 200 &&
            custom_info.deform_on &&
            custom_info.npieces == 4
        # Test warnings
        negΔt_str = "Fracturing can't occur with negative timesteps. Turning \
            fracturing off."
        nocriteria_str = "Fracturing can't occur on with NoFracture criteria. \
            Turning fracturing off."
        deform_str= "Deformation can't occur on without fracturing. Turning \
            deformation off."
        npiece_str = "Fracturing can't occur on with npieces < 2 as this won't \
            split floe. Turning fracturing off."
        # Test negative timestep provided -> can't fracture
        negΔt_info = @test_logs (:warn, negΔt_str) (:warn, deform_str) FractureSettings(
                fractures_on = true,
                criteria = criteria,
                Δt = -200,
                deform_on = true,
            )
        @test !negΔt_info.fractures_on &&
            negΔt_info.Δt == -200 &&
            !negΔt_info.deform_on
        # Test no fracture criteria provided -> can't fracture
        nocriteria_info = @test_logs (:warn, nocriteria_str) (:warn, deform_str) FractureSettings(
                fractures_on = true,
                criteria = NoFracture(),
                Δt = 200,
                deform_on = true,
            )
        @test !nocriteria_info.fractures_on && !nocriteria_info.deform_on
        # Test npieces == 1 -> won't fracture
        small_npiece_info = @test_logs (:warn, npiece_str) (:warn, deform_str) FractureSettings(
            fractures_on = true,
            criteria = criteria,
            Δt = 200,
            deform_on = true,
            npieces = 1,
        )
        @test !small_npiece_info.fractures_on && !small_npiece_info.deform_on
    end

    @testset "SimplificationSettings" begin
        # Test defaults
        default_info = SimplificationSettings()
        @test default_info.dissolve_on &&
            default_info.min_floe_area == 1e6 &&
            default_info.smooth_vertices_on &&
            default_info.max_vertices == 30 &&
            default_info.Δt_smooth == 20

        # Test warnings
        neg_min_str = "If the minimum floe area is less than or equal to 0, no floes will be dissolved. Turning dissolve floes off."
        neg_min_info =
            @test_logs (:warn, neg_min_str) SimplificationSettings(min_floe_area = 0.0)
        @test !neg_min_info.dissolve_on
        
        neg_Δt_str = "Floe smoothing can't occur on a multiple of negative timesteps. Turning floe simplification off."
        neg_Δt_info = @test_logs (:warn, neg_Δt_str) SimplificationSettings(Δt_smooth = -20)
        @test !neg_Δt_info.smooth_vertices_on
    end
end