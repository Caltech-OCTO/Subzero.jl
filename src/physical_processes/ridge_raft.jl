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
        shrinking_floe,
        growing_floe,
        vol,
        coupling_settings,
        consts,
        rng,  
    )

Removes area/volume of overlap from floe that loses area during ridging/rafting
Inputs:
    shrinking_floe      <Union{Floe, LazyRow{Floe}}> floe that loses area
    growing_floe        <Union{Floe, LazyRow{Floe}}> floe that subsumes area
    vol                 <AbstractFloat> volume of area being removed
    coupling_settings   <CouplingSettings>  simulation's coupling settings
    consts              <Constants>  simulation constants
    rng                 <AbstractRNG> random number generator
"""
function remove_floe_overlap!(
    shrinking_floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    shrinking_floeidx,
    growing_floe_coords,
    vol,
    pieces_buffer,
    coupling_settings,
    consts,
    rng,  
) where {FT <: AbstractFloat}
    # Find new floe shape
    new_floe_poly = LG.difference(
        LG.Polygon(shrinking_floe.coords),
        LG.Polygon(growing_floe_coords),
    )
    new_floe_area = LG.area(new_floe_poly)
    regions = get_polygons(new_floe_poly)
    # Update floe and create new floes from ridging
    for i in eachindex(regions)
        new_coords = find_poly_coords(regions[i])::PolyVec{FT}
        new_poly = LG.Polygon(rmholes(new_coords))
        new_mass = LG.area(new_poly) / new_floe_area *
            (shrinking_floe.mass - vol * consts.ρi)
        mass_tmp = shrinking_floe.mass
        moment_tmp = shrinking_floe.moment
        x_tmp, y_tmp = shrinking_floe.centroid
        if i == 1  # Update existing floe
            # Update floe based on new shape
            replace_floe!(
                shrinking_floe,
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
            buffer_idx = (i - 1) * shrinking_floeidx
            pieces_buffer[buffer_idx] = deepcopy_floe(shrinking_floe)
            # Update floe based on new shape
            replace_floe!(
                LazyRow(pieces_buffer, buffer_idx),
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
                LazyRow(pieces_buffer, buffer_idx),
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
    floes,
    idx1,
    idx2,
    overlap_area,
    pieces_buffer,
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
            LazyRow(floes, idx1),
            vol,
            simp_settings,
            consts,
            Δt,
        )
        extra_floes = remove_floe_overlap!(
            LazyRow(floes, idx2),
            idx2,
            floes.coords[idx1],
            vol,
            pieces_buffer,
            coupling_settings,
            consts,
            rng,  
        )
    elseif (f1_h && f2_h) ||  (!f1_h && f2_h)
        #=
        Either both floes are over max height and we randomly pick floe 2 to
        "subsumes" the extra area, or only floe 2 is over max height and it
        gets the extra area
        =#
        vol = overlap_area * floe1.height
        add_floe_volume!(
            LazyRow(floes, idx2),
            vol,
            simp_settings,
            consts,
            Δt,
        )
        extra_floes = remove_floe_overlap!(
            LazyRow(floes, idx1),
            idx1,
            floes.coords[idx2],
            vol,
            pieces_buffer,
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
    floe1_idx,
    domain_element::Union{
        CollisionBoundary,
        CompressionBoundary,
        TopographyElement,
    },
    pieces_buffer,
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
            floe1_idx,
            domain_element.coords,
            vol1,
            pieces_buffer,
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

function floe_domain_raft_ridge!(
    floe1,
    floe1_idx,
    domain_element::Union{
        PeriodicBoundary,
        OpenBoundary,
    },
    pieces_buffer,
    simp_settings,
    consts,
)
    return
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
    floes,
    idx1,
    idx2,
    overlap_area,
    pieces_buffer,
    simp_settings,
    consts,
    Δt,
)
    # Based on height ratio, pick which floe subsumes shares area
    if rand(1) >= 1/(1 + (floe1.height/floe2.height))  # Floe 1 subsumes
        vol = overlap_area * floe2.height
        # Add extra area/volume to floe 1
        add_floe_volume!(
            LazyRow(floes, idx1),
            vol,
            simp_settings,
            consts,
            Δt,
        )
        # Change shape of floe 2
        extra_floes = remove_floe_overlap!(
            LazyRow(floes, idx2),
            idx2,
            floes.coords[idx1],
            vol,
            pieces_buffer,
            coupling_settings,
            consts,
            rng,  
        )
        return extra_floes 
    else  # Floe 2 subsumes
        vol = overlap_area * floe1.height
        # Add extra area/volume to floe 2
        add_floe_volume!(
            LazyRow(floes, idx2),
            vol,
            simp_settings,
            consts,
            Δt,
        )
        # Change shape of floe 1
        extra_floes = remove_floe_overlap!(
            LazyRow(floes, idx1),
            idx1,
            floes.coords[idx2],
            vol,
            pieces_buffer,
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
    idx1,
    domain_element::AbstractDomainElement,
    pieces_buffer,
    simp_settings,
    consts,
) = floe_domain_ridge!(
    floe1,
    idx1,
    domain_element,
    pieces_buffer,
    simp_settings,
    consts,
)

function timestep_ridging!(
    floes,
    domain,
    ridgeraft_settings::RidgeRaftSettings{FT},
    pieces_buffer,
    simp_settings,
    consts,
    Δt,
) where {FT <: AbstractFloat}
    pieces_buffer = Vector{Union{nothing, PolyVec{FT}}}(nothing, length(floes) * 2)
    for i in eachindex(floes)
        #=
            Floe is active in sim, hasn't broken, interacted with other floes,
            isn't too thick, random check meets probability to ridge
        =#
        if (
            floes.status[i].tag == active &&
            isnothing(pieces_buffer[i]) &&
            !isempty(floes.interactions[i]) &&
            floes.height[i] <= ridgeraft_settings.max_floe_ridge_height &&
            rand() <= ridgeraft_settings.ridge_probability
        )        
            # TODO: Need to make sure we don't have repeats 
            for jrow in axes(floes.interactions[i], 1)
                j = Int(floes.interactions[jrow, floeidx])
                #=
                    Ridge between two floes (j > 0) if neither are broken and if
                    both floes aren't too thick
                =#
                if (
                    j > 0 &&
                    isnothing(pieces_buffer[i]) &&
                    isnothing(pieces_buffer[j]) &&
                    floes.height[i] <= ridgeraft_settings.max_floe_ridge_height &&
                    floes.height[j] <= ridgeraft_settings.max_floe_ridge_height
                )
                    floe_floe_ridge!(
                        floes,
                        i,
                        j,
                        floe.interactions[i, overlap],
                        pieces_buffer,
                        ridgeraft_settings,
                        simp_settings,
                        consts,
                        Δt,
                    )
                #=
                    Ridge between floe and domain element (j < 0) if floe isn't
                    broken and if the floe isn't too thick
                =#
                elseif (
                    j < 0 &&
                    isnothing(pieces_buffer[i]) &&
                    floes.height[i] <= ridgeraft_settings.max_domain_ridge_height
                )
                    floe_domain_ridge!(
                        LazyRow(floes, i),
                        i,
                        get_domain_element(domain, j),
                        pieces_buffer,
                        simp_settings,
                        consts,
                    )
                end
            end
        end
    end
    return
end

function timestep_rafting!(
    floes,
    domain,
    ridgeraft_settings::RidgeRaftSettings{FT},
    pieces_buffer,
    simp_settings,
    consts,
    Δt,
) where {FT <: AbstractFloat}
    for i in eachindex(floes)
        #=
            Floe is active in sim, hasn't broken, interacted with other floes,
            isn't too thick, random check meets probability to raft
        =#
        if (
            floes.status[i].tag == active &&
            isnothing(pieces_buffer[i]) &&
            !isempty(floes.interactions[i]) &&
            floes.height[i] <= ridgeraft_settings.max_floe_raft_height &&
            rand() <= ridgeraft_settings.raft_probability
        )
            for jrow in axes(floes.interactions[i], 1)
                j = Int(floes.interactions[jrow, floeidx])
                #=
                    Raft between two floes (j > 0) if neither are broken and if
                    both floes aren't too thick
                =#
                if (
                    j > 0 &&
                    isnothing(pieces_buffer[i]) &&
                    isnothing(pieces_buffer[j]) &&
                    floes.height[i] <= ridgeraft_settings.max_floe_raft_height &&
                    floes.height[j] <= ridgeraft_settings.max_floe_raft_height
                )
                    floe_floe_ridge!(
                        floes,
                        i,
                        j,
                        floe.interactions[i, overlap],
                        pieces_buffer,
                        ridgeraft_settings,
                        simp_settings,
                        consts,
                        Δt,
                    )
                #=
                    Raft between floe and domain element (j < 0) if floe isn't
                    broken and if the floe isn't too thick
                =#
                elseif (
                    j < 0 &&
                    isnothing(pieces_buffer[i]) &&
                    floes.height[i] <= ridgeraft_settings.max_domain_raft_height
                )
                    floe_domain_raft!(
                        LazyRow(floes, i),
                        i,
                        get_domain_element(domain, j),
                        pieces_buffer,
                        simp_settings,
                        consts,
                    )
                end

            end
        end
    end

    return
end
