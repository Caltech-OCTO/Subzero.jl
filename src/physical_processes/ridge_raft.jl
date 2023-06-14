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
)
    init_height = floe.height
    floe.height += vol/floe.area
    if floe.height > simp_settings.max_floe_height
        floe.height = simp_settings.max_floe_height
    end
    floe.mass += vol * consts.ρi
    floe.moment *= floe.height/init_height
end

function floe_floe_ridge!(
    floe1,
    floe2,
    overlap_area,
    ridgeraft_settings,
    simp_settings,
    consts,
)
    vol1 = overlap_area * floe1.height
    vol2 = overlap_area * floe2.height

    if floe1.height >= ridgeraft_settings.hc && floe2.height >= ridgeraft_settings.hc
        if rand(1) >= 1/(1 + (floe1.height/floe2.height))
            add_floe_volume!(
                floe1,
                overlap_area * floe2.height,
                simp_settings,
                consts,
            )
             # ridge value switch 
        else
            add_floe_volume!(
                floe2,
                overlap_area * floe1.height,
                simp_settings,
                consts,
            )
             # ridge value switch 
        end
    elseif floe1.height >= ridgeraft_settings.hc && floe2.height < ridgeraft_settings.hc
        add_floe_volume!(
            floe1,
            overlap_area * floe2.height,
            simp_settings,
            consts,
        )
         # ridge value switch 
    elseif floe1.height < ridgeraft_settings.hc && floe2.height >= ridgeraft_settings.hc
        add_floe_volume!(
            floe2,
            overlap_area * floe1.height,
            simp_settings,
            consts,
        )
         # ridge value switch 
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