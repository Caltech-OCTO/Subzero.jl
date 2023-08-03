"""
Functions needed for ridging and rafting between floes and boundaries. 
"""

"""
    add_floe_volume!(
        floe,
        vol,
        simp_settings,
        consts,
        Δt,
    )
Add volume to existing floe and update fields that depend on the volume.
Inputs:
    floe            <Union{Floe, LazyRow{Floe}}> floe to add volume to
    vol             <AbstractFloat> volume to add to floe
    simp_settings   <SimplificationSettings> simulation simplification settings
    consts          <Constants>  simulation's constants
    Δt              <Int> length of timestep in seconds
Outputs:
    Nothing. Floe's fields are updated to reflect increase in volume.
"""
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

"""
    remove_floe_overlap!(
        keepfloe,
        floediff,
        vol,
        coupling_settings,
        consts,
        rng,  
    )

Removes area/volume of overlap from floe that loses area during ridging/rafting
Inputs:
    keepfloe            <Union{Floe, LazyRow{Floe}}> floe that loses area
    floediff            <Union{Floe, LazyRow{Floe}}> floe that subsumes area
    vol                 <AbstractFloat> volume of area being removed
    coupling_settings   <CouplingSettings>  simulation's coupling settings
    consts              <Constants>  simulation constants
    rng                 <AbstractRNG> random number generator
"""
function remove_floe_overlap!(
    keepfloe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    floediff,
    vol,
    coupling_settings,
    consts,
    rng,  
) where {FT <: AbstractFloat}
    # Find new floe shape
    new_floe_poly = LG.difference(
        LG.Polygon(keepfloe.coords),
        LG.Polygon(floediff.coords),
    )
    new_floe_area = LG.area(new_floe_poly)
    regions = get_polygons(new_floe_poly)
    new_floes = StructVector{Floe}(undef)
    n_new = 0
    # Update floe and create new floes from ridging
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
        else  # Create new floes if floe was split into multiple pieces
            push!(new_floes, deepcopy_floe(keepfloe))
            n_new += 1
            # Update floe based on new shape
            replace_floe!(
                LazyRow(new_floes, n_new),
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
                LazyRow(new_floes, n_new),
            )
        end
    end
end

"""
    floe_floe_ridge!(
        floe1,
        floe2,
        overlap_area,
        ridgeraft_settings,
        simp_settings,
        consts,
        Δt,
    )
Ridge two floes, updating both in-place and returning any new floes that
resulting from the ridging event.
Inputs:
    floe1               <Union{Floe, LazyRow{Floe}}> floe 1 involved in ridging
    floe2               <Union{Floe, LazyRow{Floe}}> floe 2 involved in ridging
    overlap_area        <AbstractFloat> overlap area between floes
    ridgeraft_settings  <RidgeRaftSettings> simulation's settings for ridging
                            and rafting
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    consts              <Constants> simulation's constants
    Δt                  <Int> simulation timestep in seconds
Outputs:
    Updates floe1 and floe2 and returns any new floes created by ridging
"""
function floe_floe_ridge!(
    floe1,
    floe2,
    overlap_area,
    ridgeraft_settings,
    simp_settings,
    consts,
    Δt,
)
    # Heights of floes determine which floe subsumes shared area
    f1_h = floe1.height >= ridgeraft_settings.min_ridge_height
    f2_h = floe2.height >= ridgeraft_settings.min_ridge_height

    extra_floes = if(
        (f1_h && f2_h && rand() >= 1/(1 + (floe1.height/floe2.height))) ||
        (f1_h && !f2_h)
    )
        #=
        Either both floes are over max height and we randomly pick floe 1 to
        "subsume" the extra area, or only floe 1 is over max height and it
        gets the extra area
        =#
        vol = overlap_area * floe2.height
        add_floe_volume!(
            floe1,
            vol,
            simp_settings,
            consts,
            Δt,
        )
        extra_floes = remove_floe_overlap!(
            floe2,
            floe1,
            vol,
            coupling_settings,
            consts,
            rng,  
        )
        return extra_floes
    elseif (f1_h && f2_h) ||  (!f1_h && f2_h)
        #=
        Either both floes are over max height and we randomly pick floe 2 to
        "subsumes" the extra area, or only floe 2 is over max height and it
        gets the extra area
        =#
        vol = overlap_area * floe1.height
        add_floe_volume!(
            floe2,
            vol,
            simp_settings,
            consts,
            Δt,
        )
        extra_floes = remove_floe_overlap!(
            floe1,
            floe2,
            vol,
            coupling_settings,
            consts,
            rng,  
        )
        return extra_floes
    end
    return nothing
end

"""
    floe_domain_ridge!(
        floe1,
        domain_element,
        simp_settings,
        consts,
    )

Ridge a floe against a boundary or a topography element and return any excess
floes created by the ridging.
Inputs:
    floe1           <Union{Floe, LazyRow{Floe}}> floe ridging
    domain_element  <AbstractDomainElement> boundary or topography element
    simp_settings   <SimplificationSettings> simulation simplification settings
    consts          <Constants> simulation constants
Outputs:
    floe1 is updated with new shape. If any new floes are created by ridging
    they are returned, else nothing.
"""
function floe_domain_ridge!(
    floe1::Union{Floe, LazyRow{Floe}},
    domain_element::AbstractDomainElement,
    simp_settings,
    consts,
)
    # Calculate new floe shape and area/volume lost
    ridged_floe = LG.difference(
        LG.Polygon(floe1.coords),
        LG.Polygon(domain_element.coords)    
    )
    ridged_area = LG.area(ridged_floe)
    vol1 = ridged_area * floe1.height
    # Ridge floe if area is large enough, else reduce area so it will be removed
    if ridged_area > simp_settings.min_floe_area  # boundary subsumes overlap
        extra_floes = remove_floe_overlap!(
            floe1,
            domain_element,
            vol1,
            coupling_settings,
            consts,
            rng,  
        )
        return extra_floes
    else  # remaining floe is too small, decrease area/mass to remove
        floe1.area = ridged_area
        floe1.mass = ridged_area * floe1.height * consts.ρi
        # do I need to do more so that this doesn't get used in other sections?
        # maybe we do need dissolved tag...
        return nothing
    end
    return nothing
end

"""
    floe_floe_raft!(
        floe1::F,
        floe2::F,
        overlap_area,
        simp_settings,
        consts,
        Δt,
    )

Raft two floes, updating both in-place and returning any new floes that
resulting from the rafting event.
Inputs:
    floe1               <Union{Floe, LazyRow{Floe}}> floe 1 involved in rafting
    floe2               <Union{Floe, LazyRow{Floe}}> floe 2 involved in rafting
    overlap_area        <AbstractFloat> overlap area between floes
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    consts              <Constants> simulation's constants
    Δt                  <Int> simulation timestep in seconds
Outputs:
    Updates floe1 and floe2 and returns any new floes created by rafting
"""
function floe_floe_raft!(
    floe1::Union{Floe, LazyRow{Floe}},
    floe2::Union{Floe, LazyRow{Floe}},
    overlap_area,
    simp_settings,
    consts,
    Δt,
)
    # Based on height ratio, pick which floe subsumes shares area
    if rand(1) >= 1/(1 + (floe1.height/floe2.height))  # Floe 1 subsumes
        vol = overlap_area * floe2.height
        # Add extra area/volume to floe 1
        add_floe_volume!(
            floe1,
            vol,
            simp_settings,
            consts,
            Δt,
        )
        # Change shape of floe 2
        extra_floes = remove_floe_overlap!(
            floe2,
            floe1,
            vol,
            coupling_settings,
            consts,
            rng,  
        )
        return extra_floes 
    else  # Floe 2 subsumes
        vol = overlap_area * floe1.height
        # Add extra area/volume to floe 2
        add_floe_volume!(
            floe2,
            vol,
            simp_settings,
            consts,
            Δt,
        )
        # Change shape of floe 1
        extra_floes = remove_floe_overlap!(
            floe1,
            floe2,
            vol,
            coupling_settings,
            consts,
            rng,  
        )
        return extra_floes
    end
    return nothing
end

"""
    floe_domain_raft(
        floe1::Union{Floe, LazyRow{Floe}},
        domain_element::AbstractDomainElement,
        simp_settings,
        consts,
    )

Raft a floe against a boundary or a topography element and return any excess
floes created by the rafting. This is equivalent to ridging.
Inputs:
    floe1           <Union{Floe, LazyRow{Floe}}> floe rafting
    domain_element  <AbstractDomainElement> boundary or topography element
    simp_settings   <SimplificationSettings> simulation simplification settings
    consts          <Constants> simulation constants
Outputs:
    floe1 is updated with new shape. If any new floes are created by rafting
    they are returned, else nothing.
"""
floe_domain_raft(
    floe1::Union{Floe, LazyRow{Floe}},
    domain_element::AbstractDomainElement,
    simp_settings,
    consts,
) = floe_domain_ridge!(
    floe1,
    domain_element,
    simp_settings,
    consts,
)

function timestep_ridging_rafting!(
    floes,
    domain,
    ridgeraft_settings::RidgeRaftSettings{FT},
    pieces_buffer,
    simp_settings,
    consts,
    Δt,
) where {FT <: AbstractFloat}
    for i in eachindex(floes)
        # ridging and rafting only happens with floes that collided
        if (
            floes.status[i].tag != remove &&
            !isempty(floes.interactions[i])
        )        
            # Ridge floe if it meets ridge probabilities and maximum height
            if (floes.height[i] <= ridgeraft_settings.max_ridge_height &&
                rand() <= ridgeraft_settings.ridge_probability
            )
                for j in axes(floes.interactions[i], 1)
                    # Ridge between two floes


                    # Ridge between floe and domain element --> make sure don't do for periodic
                end
            end
            # Raft floe if it meets raft possibility and maximum height
            if (floes.height[i] <= ridgeraft_settings.max_raft_height &&
                rand() <= ridgeraft_settings.raft_probability
            )
                # Raft between two floes

                # Raft between floe and domain element --> make sure don't do for periodic

            end
        end
    end

end
