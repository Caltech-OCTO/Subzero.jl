@testset "ProcessInfo" begin
    @testset "CouplingInfo" begin
        # Test default
        default_info = CouplingInfo()
        @test default_info.coupling_on && default_info.Δt == 10 && !default_info.calc_ocnτ
        # Test custom correct
        custom_info = CouplingInfo(false, 20, false)
        @test !custom_info.coupling_on && custom_info.Δt == 20 && !custom_info.calc_ocnτ
        # Test negative timestep -> turns coupling off so turn calc_ocnτ off 
        warn_str1 = "Coupling can't occur on a multiple of negative timesteps. Turning coupling off."
        warn_str2 = "Can't calculate stresses on ocean from ice and atmosphere without coupling. Turning calc_ocnτ off."
        neg_Δt_info = @test_logs (:warn, warn_str1) (:warn, warn_str2) CouplingInfo(true, -20, true)
        @test !neg_Δt_info.coupling_on && neg_Δt_info.Δt == -20 && !neg_Δt_info.calc_ocnτ
    end

    @testset "CollisionInfo" begin
        # Test default
        default_info = CollisionInfo()
        @test default_info.collisions_on && default_info.floe_floe_max_overlap == 0.55 && default_info.floe_domain_max_overlap == 0.75
        # Test custom false
        custom_info = CollisionInfo(false, -0.1, 2.0)
        @test !custom_info.collisions_on && custom_info.floe_floe_max_overlap == -0.1 && custom_info.floe_domain_max_overlap == 2
        # Test overlaps out of bounds
        floe2large_str = "The maximum collisin overlap between floes can't be greater than 1. Setting to 1."
        floe2small_str = "The maximum collisin overlap between floes can't be less than 0. Setting to 0."
        domain2large_str = "The maximum collisin overlap between floes and the domain can't be greater than 1. Setting to 1."
        domain2small_str = "The maximum collisin overlap between floes and the domain can't be less than 0. Setting to 0."
        out_of_bounds = @test_logs (:warn, floe2large_str) (:warn, domain2small_str) CollisionInfo(true, 1.2, -0.1)
        @test out_of_bounds.collisions_on && out_of_bounds.floe_floe_max_overlap == 1.0 && out_of_bounds.floe_domain_max_overlap == 0.0
        
        out_of_bounds = @test_logs (:warn, floe2small_str) (:warn, domain2large_str) CollisionInfo(true, -0.1, 1.2)
        @test out_of_bounds.collisions_on && out_of_bounds.floe_floe_max_overlap == 0.0 && out_of_bounds.floe_domain_max_overlap == 1.0

    end
    @testset "FractureInfo" begin
        default_info = FractureInfo()
        @test !default_info.fractures_on && default_info.criteria isa NoFracture && default_info.Δt == 0 && !default_info.deform_on
        criteria = HiblerYieldCurve(0.0, 0.0, [[[0.0]]])
        custom_info = FractureInfo(true, criteria, 200, true)
        @test custom_info.fractures_on && custom_info.criteria isa HiblerYieldCurve && custom_info.Δt == 200 && custom_info.deform_on
        negΔt_str = "Fracturing can't occur on a multiple of negative timesteps. Turning fracturing off."
        nocriteria_str = "Fracturing can't occur on with NoFracture criteria. Turning fracturing off."
        deform_str= "Deformation can't occur on without fracturing. Turning deformation off."
        negΔt_info = @test_logs (:warn, negΔt_str) (:warn, deform_str)  FractureInfo(true, criteria, -200, true)
        @test !negΔt_info.fractures_on && negΔt_info.Δt == -200 && !negΔt_info.deform_on
        nocriteria_info = @test_logs (:warn, nocriteria_str) (:warn, deform_str) FractureInfo(true, NoFracture(), 200, true)
        @test !nocriteria_info.fractures_on && !negΔt_info.deform_on
    end
end