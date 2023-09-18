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
    floes,
    idx,
    vol,
    simp_settings,
    consts,
    Δt,
)
    # Find parent floe (could be current floe, or current floe could be ghost)
    if floes.ghost_id[idx] != 0
        idx = findfirst(x -> x == floes.id[idx], floes.id)
    end
    # Update floe height, mass, and moment of intertia due to volume change
    init_height = floes.height[idx]
    floes.height[idx] += vol/floes.area[idx]
    if floes.height[idx] > simp_settings.max_floe_height
        floes.height[idx] = simp_settings.max_floe_height
    end
    floes.mass[idx] += vol * consts.ρi
    floes.moment[idx] *= floes.height[idx] / init_height
    # Update all ghost floes of parent
    for gidx in floes.ghosts[idx]
        floes.height[gidx] = floes.height[idx]
        floes.mass[gidx] = floes.mass[idx]
        floes.moment[gidx] = floes.moment[idx]
    end
    # TODO: Update velocities and accelerations to conserve momentum
    return
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
    floes::StructArray{<:Floe{FT}},
    shrinking_idx,
    growing_floe_coords,
    pieces_buffer,
    max_floe_id,
    broken,
    coupling_settings,
    consts,
    rng,  
) where {FT <: AbstractFloat}
    # Find new floe shape and regions
    new_floe_poly = LG.difference(
        LG.Polygon(floes.coords[shrinking_idx]),
        LG.Polygon(growing_floe_coords),
    )
    regions = get_polygons(new_floe_poly)
    nregions = length(regions)
    # Calculate changes in mass / area
    new_floe_area = LG.area(new_floe_poly)
    transfer_vol = (floes.area[shrinking_idx] - new_floe_area) * floes.height[shrinking_idx]
    new_floe_mass =  (floes.mass[shrinking_idx] - transfer_vol * consts.ρi)
    # Reset shrinking index to parent floe and determine floe shift
    parent_Δx = FT(0)
    parent_Δy = FT(0)
    if floes.ghost_id[shrinking_idx] != 0
        parent_idx = findfirst(x -> x == floes.id[shrinking_idx], floes.id)
        parent_Δx = floes.centroid[parent_idx][1] - floes.centroid[shrinking_idx][1]
        parent_Δy = floes.centroid[parent_idx][2] - floes.centroid[shrinking_idx][2]
        shrinking_idx = parent_idx
    end
    parent_centroid = floes.centroid[shrinking_idx]
    # Update existing floes/ghosts regions
    for i in 1:nregions
        new_coords = find_poly_coords(regions[i])::PolyVec{FT}
        translate!(  # shift region coords to parent floe location
            new_coords,
            parent_Δx,
            parent_Δy,
        )
        rmholes!(new_coords)  # remove holes in floe
        new_poly = LG.Polygon(new_coords)  # parent floes new region polygon
        new_mass = (LG.area(new_poly) / new_floe_area) * new_floe_mass
        buffer_length = length(pieces_buffer)
        if i == 1
            replace_floe!(  # replace parent floe
                LazyRow(floes, shrinking_idx),
                new_poly,
                new_mass,
                coupling_settings,
                consts,
                rng,
            )
            if nregions == 1
                # if floe doesn't break, propograte changes to ghosts
                for gidx in floes.ghosts[shrinking_idx]  # update any ghosts
                    # find shift to ghost floe location
                    g_Δx = floes.centroid[gidx][1] - parent_centroid[1]
                    g_Δy = floes.centroid[gidx][2] - parent_centroid[2]
                    # replace ghost floe
                    replace_floe!(
                        LazyRow(floes, gidx),
                        new_poly,
                        new_mass,
                        coupling_settings,
                        consts,
                        rng,
                    )
                    # shift ghost floe
                    translate!(floes.coords[gidx], g_Δx, g_Δy)
                    floes.centroid[gidx][1] += g_Δx
                    floes.centroid[gidx][1] += g_Δy
                end
            else
                # if floe breaks, mark floe and ghosts as broken
                broken[shrinking_idx] = true
                for gidx in floes.ghosts[shrinking_idx]
                    broken[gidx] = true
                    floes.status[gidx].tag = remove
                end
                # Update floe identifiers
                empty!(floes.ghosts[shrinking_idx])
                push!(floes.parent_ids[shrinking_idx], floes.id[shrinking_idx])
                max_floe_id += 1
                floes.id[shrinking_idx] = max_floe_id
            end
        else  # >1 region, so floe must break and add pieces to buffer
            push!(
                pieces_buffer,
                deepcopy_floe(LazyRow(floes, shrinking_idx))
            )
            buffer_length += 1
            replace_floe!(
                LazyRow(pieces_buffer, buffer_length),
                new_poly,
                new_mass,
                coupling_settings,
                consts,
                rng,
            )
            max_floe_id += 1
            pieces_buffer.id[buffer_length] = max_floe_id
        end
    end
    return transfer_vol, max_floe_id
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
    floes::StructArray{<:Floe{FT}},
    idx1,
    idx2,
    pieces_buffer,
    max_floe_id,
    broken,
    ridgeraft_settings,
    coupling_settings,
    simp_settings,
    consts,
    Δt,
    rng,
) where {FT}
    # Heights of floes determine which floe subsumes shared area
    f1_h = floes.height[idx1] >= ridgeraft_settings.min_ridge_height
    f2_h = floes.height[idx2] >= ridgeraft_settings.min_ridge_height
    vol = FT(0)
    if(
        (f1_h && f2_h && rand() >= 1/(1 + (floes.height[idx1]/floes.height[idx2]))) ||
        (f1_h && !f2_h)
    )
        #=
        Either both floes are over min height and we randomly pick floe 1 to
        "subsume" the extra area, or only floe 1 is over min height and it
        gets the extra area
        =#
        vol, max_floe_id = remove_floe_overlap!(
            floes,
            idx2,
            floes.coords[idx1],
            pieces_buffer,
            max_floe_id,
            broken,
            coupling_settings,
            consts,
            rng,  
        )
        add_floe_volume!(
            floes,
            idx1,
            vol,
            simp_settings,
            consts,
            Δt,
        )
    elseif (f1_h && f2_h) ||  (!f1_h && f2_h)
        #=
        Either both floes are over max height and we randomly pick floe 2 to
        "subsumes" the extra area, or only floe 2 is over max height and it
        gets the extra area
        =#
        vol, max_floe_id = remove_floe_overlap!(
            floes,
            idx1,
            floes.coords[idx2],
            pieces_buffer,
            max_floe_id,
            broken,
            coupling_settings,
            consts,
            rng,  
        )
        add_floe_volume!(
            floes,
            idx2,
            vol,
            simp_settings,
            consts,
            Δt,
        )
    end
    return max_floe_id
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
    floes,
    idx,
    domain_element::Union{<:AbstractDomainElement, LazyRow{<:AbstractDomainElement}},
    pieces_buffer,
    max_floe_id,
    broken,
    coupling_settings,
    simp_settings,
    consts,
    rng,
)
    _, max_floe_id = remove_floe_overlap!(
        floes,
        idx,
        domain_element.coords,
        pieces_buffer,
        max_floe_id,
        broken,
        coupling_settings,
        consts,
        rng,  
    )
    return max_floe_id
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
    floes::StructArray{<:Floe{FT}},
    idx1,
    idx2,
    pieces_buffer,
    max_floe_id,
    broken,
    coupling_settings,
    simp_settings,
    consts,
    Δt,
    rng,
) where {FT}
    # Based on height ratio, pick which floe subsumes shares area
    vol = FT(0)
    if rand(rng) >= 1/(1 + (floes.height[idx1]/floes.height[idx2]))  # Floe 1 subsumes
        # Change shape of floe 2
        vol, max_floe_id = remove_floe_overlap!(
            floes,
            idx2,
            floes.coords[idx1],
            pieces_buffer,
            max_floe_id,
            broken,
            coupling_settings,
            consts,
            rng,  
        )
        # Add extra area/volume to floe 1
        add_floe_volume!(
            floes,
            idx1,
            vol,
            simp_settings,
            consts,
            Δt,
        )
    else  # Floe 2 subsumes
        # Change shape of floe 1
        vol, max_floe_id = remove_floe_overlap!(
            floes,
            idx1,
            floes.coords[idx2],
            pieces_buffer,
            max_floe_id,
            broken,
            coupling_settings,
            consts,
            rng,  
        )
        # Add extra area/volume to floe 2
        add_floe_volume!(
            floes,
            idx2,
            vol,
            simp_settings,
            consts,
            Δt,
        )
    end
    return max_floe_id
end

"""
    floe_domain_raft!(
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
floe_domain_raft!(
    floes,
    idx1,
    domain_element::Union{AbstractDomainElement, LazyRow{<:AbstractDomainElement}},
    pieces_buffer,
    max_floe_id,
    broken,
    coupling_settings,
    simp_settings,
    consts,
    rng,
) = floe_domain_ridge!(
    floes,
    idx1,
    domain_element,
    pieces_buffer,
    max_floe_id,
    broken,
    coupling_settings,
    simp_settings,
    consts,
    rng,
)

"""
    timestep_ridging_rafting!(
        floes,
        domain,
        ridgeraft_settings::RidgeRaftSettings{FT},
        pieces_buffer,
        simp_settings,
        consts,
        Δt,
    )

Ridge and raft floes that meet probability and height criteria.
Inputs:
    floes               <StructArray{Floe}> simulation's list of floes
    n_init_floes        <Int> number of floes prior to adding ghosts
    domain              <Domain> simulation's domain
    ridgeraft_settings  <RidgeRaftSettings> ridge/raft settings
    pieces_buffer       <Vector{Nothing, ???}> buffer to hold extra floe pieces
    simp_settings       <SimplificationSettings> simplification settings
    consts              <Consts> simulation's constants
    Δt                  <Int> length of timestep in seconds
Outputs:
    Updates floes post ridging and rafting and adds any new pieces to the pieces
    buffer to be made into new floes.
"""
function timestep_ridging_rafting!(
    floes,
    pieces_list,
    n_init_floes,
    domain,
    max_floe_id,
    ridgeraft_settings::RidgeRaftSettings{FT},
    coupling_settings,
    simp_settings,
    consts,
    Δt,
    rng = Xoshiro(),
) where {FT <: AbstractFloat}
    broken = fill(false, length(floes))
    interactions_list = zeros(Int, size(floes.interactions[1], 1))
    for i in 1:n_init_floes
        #=
            Floe is active in sim, hasn't broken, interacted with other floes,
            isn't too thick, and meets probability check to ridge or raft
        =#
        ridge = floes.height[i] <= ridgeraft_settings.max_floe_ridge_height &&
            rand() <= ridgeraft_settings.ridge_probability
        raft = floes.height[i] <= ridgeraft_settings.max_floe_raft_height &&
            rand() <= ridgeraft_settings.raft_probability
        if (
            (ridge || raft) &&
            floes.status[i].tag == active &&
            !broken[i] &&
            floes.num_inters[i] > 0
        )
            ninters = 0
            # Find unique interactions that haven't already been calculated
            for row in 1:floes.num_inters[i]
                j = Int(floes.interactions[i][row, floeidx])
                min_area = j > 0 ?
                    min(floes.area[i], floes.area[j]) :  # floe-floe interaction
                    floes.area[i] # floe-domain interaction
                if (
                    ((j > i  && !broken[j]) || j < 0) &&
                    1e-6 < floes.interactions[i][row, overlap]/min_area < 0.95 &&
                    !(j in interactions_list)
                )
                    ninters += 1
                    if ninters <= length(interactions_list)
                        interactions_list[ninters] = j
                    else
                        push!(interactions_list, j)
                    end
                end     
            end 
            for k in 1:ninters
                # if floe-floe interaction
                j = interactions_list[k]
                if j > 0 && !broken[i] && !broken[j]
                    if (  # if both floes aren't too thick -> ridge
                        ridge &&
                        floes.height[i] <= ridgeraft_settings.max_floe_ridge_height &&
                        floes.height[j] <= ridgeraft_settings.max_floe_ridge_height
                    )
                        max_floe_id = floe_floe_ridge!(
                            floes,
                            i,
                            j,
                            pieces_list,
                            max_floe_id,
                            broken,
                            ridgeraft_settings,
                            coupling_settings,
                            simp_settings,
                            consts,
                            Δt,
                            rng,
                        )
                    elseif ( # if both floes aren't too thick -> raft
                        raft &&
                        floes.height[i] <= ridgeraft_settings.max_floe_raft_height &&
                        floes.height[j] <= ridgeraft_settings.max_floe_raft_height
                    )
                        max_floe_id = floe_floe_raft!(
                            floes,
                            i,
                            j,
                            pieces_list,
                            max_floe_id,
                            broken,
                            coupling_settings,
                            simp_settings,
                            consts,
                            Δt,
                            rng,
                        )
                    end
                elseif j < 0 && !broken[i]  # if floe-domain interaction
                    if (  # if floe isn't too thick -> ridge
                        ridge &&
                        floes.height[i] <= ridgeraft_settings.max_domain_ridge_height
                    )
                        max_floe_id = floe_domain_ridge!(
                            floes,
                            i,
                            get_domain_element(domain, j),
                            pieces_list,
                            max_floe_id,
                            broken,
                            coupling_settings,
                            simp_settings,
                            consts,
                            rng,
                        )

                    elseif (  # if floe isn't too thick -> raft
                        raft &&
                        floes.height[i] <= ridgeraft_settings.max_domain_raft_height
                    )
                        max_floe_id = floe_domain_raft!(
                            floes,
                            i,
                            get_domain_element(domain, j),
                            pieces_list,
                            max_floe_id,
                            broken,
                            coupling_settings,
                            simp_settings,
                            consts,
                            rng,
                        )
                    end
                end
            end
        end
    end
    return max_floe_id
end
