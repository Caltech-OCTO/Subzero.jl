"""
Functions needed for ridging and rafting between floes and boundaries. 
"""

# abound == "a boundary"
# Only let one floe ridge with another floe...

function add_floe_volume!(
    floe,
    vol,
    simp_settings,
    consts,
    Δt,
)
    # Update floe height, mass, and moment of intertia due to volume change
    init_height = floe.height
    floe.height += vol/floe.area
    if floe.height > simp_settings.max_floe_height
        floe.height = simp_settings.max_floe_height
    end
    mass_tmp = floe.mass
    floe.mass += vol * consts.ρi
    moment_tmp = floe.moment
    floe.moment *= floe.height/init_height
    # Update velocities and accelerations to conserve momentum
    conserve_momentum_combination!(
        mass_tmp,
        moment_tmp,
        floe.centroid[1],
        floe.centroid[2],
        Δt,
        floe,
    )
end

function remove_floe_overlap!(
    keepfloe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    floediff,
    vol,
    coupling_settings,
    consts,
    rng,  
) where {FT <: AbstractFloat}
    new_floe_poly = LG.difference(
        LG.Polygon(keepfloe.coords),
        LG.Polygon(floediff.coords),
    )
    new_floe_area = LG.area(new_floe_poly)
    regions = get_polygons(new_floe_poly)
    new_floes = StructVector{Floe}(undef, length(regions) - 1)
    for i in eachindex(regions)
        new_coords = find_poly_coords(regions[i])::PolyVec{FT}
        new_poly = LG.Polygon(rmholes(new_coords))
        new_area = LG.area(new_poly)
        new_mass = new_area/new_floe_area * (keepfloe.mass - vol * consts.ρi)
        mass_tmp = keepfloe.mass
        moment_tmp = keepfloe.moment
        x_tmp, y_tmp = keepfloe.centroid
        if i == 1  # Update existing floe
            # Update floe based on new shape
            replace_floe!(
                keepfloe,
                new_poly,
                new_mass,
                consts,
                coupling_settings.mc_n,
                rng,
            )
            # Update velocities and accelerations to conserve momentum
            conserve_momentum_combination!(
                mass_tmp,
                moment_tmp,
                x_tmp,
                y_tmp,
                Δt,
                keep_floe,
            )
        else  # Create new floes
            new_floes[i] = deepcopy(keepfloe)
            # Update floe based on new shape
            replace_floe!(
                LazyRow(new_floes, i),
                new_poly,
                new_mass,
                consts,
                coupling_settings.mc_n,
                rng,
            )
            # Update velocities and accelerations to conserve momentum
            conserve_momentum_combination!(
                mass_tmp,
                moment_tmp,
                x_tmp,
                y_tmp,
                Δt,
                LazyRow(new_floes, i),
            )
        end
    end
end

function floe_floe_ridge!(
    floe1,
    floe2,
    overlap_area,
    ridgeraft_settings,
    simp_settings,
    consts,
    Δt,
)
    vol1 = overlap_area * floe1.height
    vol2 = overlap_area * floe2.height

    if floe1.height >= ridgeraft_settings.hc && floe2.height >= ridgeraft_settings.hc
        if rand(1) >= 1/(1 + (floe1.height/floe2.height))
            add_floe_volume!(
                floe1,
                vol2,
                simp_settings,
                consts,
                Δt,
            )
            remove_floe_overlap!(
                floe2,
                floe1,
                vol2,
                coupling_settings,
                consts,
                rng,  
            ) 
        else
            add_floe_volume!(
                floe2,
                vol1,
                simp_settings,
                consts,
                Δt,
            )
            remove_floe_overlap!(
                floe1,
                floe2,
                vol1,
                coupling_settings,
                consts,
                rng,  
            )  
        end
    elseif floe1.height >= ridgeraft_settings.hc && floe2.height < ridgeraft_settings.hc
        add_floe_volume!(
            floe1,
            vol2,
            simp_settings,
            consts,
            Δt,
        )
        remove_floe_overlap!(
            floe2,
            floe1,
            vol2,
            coupling_settings,
            consts,
            rng,  
        ) 
    elseif floe1.height < ridgeraft_settings.hc && floe2.height >= ridgeraft_settings.hc
        add_floe_volume!(
            floe2,
            vol1,
            simp_settings,
            consts,
            Δt,
        )
        remove_floe_overlap!(
            floe1,
            floe2,
            vol1,
            coupling_settings,
            consts,
            rng,  
        ) 
    end
end

function floe_domain_ridge!(
    floe1::F,
    domain_element::B,
    simp_settings,
    consts,
) where {F <: Union{Floe, LazyRow{Floe}}, B <: AbstractDomainElement}
    ridged_floe = LG.difference(
        LG.Polygon(floe1.coords),
        LG.Polygon(domain_element.coords)    
    )
    ridged_area = LG.area(ridged_floe)
    if ridged_area > simp_settings.min_floe_area
        # Add ridged floe/floes to floe list
    else
        floe1.area = ridged_area
        floe1.mass = ridged_area * floe1.height * consts.ρi
        # do I need to do more so that this doesn't get used in other sections?
        # maybe we do need dissolved tag...
    end

    return
end

function floe_floe_raft(
    floe1::F,
    floe2::F,
    overlap_area,
    ridgeraft_settings,
    simp_settings,
    consts,
) where {F <: Union{Floe, LazyRow{Floe}}}
    vol1 = overlap_area * floe1.height
    vol2 = overlap_area * floe2.height

    if rand(1) >= 1/(1 + (floe1.height/floe2.height))
        # ridge value switch 
    else
         # ridge value switch 
    end

end

floe_domain_raft(
    floe1,
    domain_element,
    simp_settings,
    consts,
) = ridge_floe(
    floe1,
    domain_element,
    simp_settings,
    consts,
)