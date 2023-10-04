"""
Functions needed for ridging and rafting between floes and boundaries. 
"""

"""
    add_floe_volume!(
        floes,
        idx,
        vol,
        simp_settings,
        consts,
        Δt,
    )
Add volume to existing floe and update fields that depend on the volume.
Inputs:
    floes           <StructArray{Frloe}> list of floes
    idx             <Int> index of floe to add volume to 
    vol             <AbstractFloat> volume to add to floe
    simp_settings   <SimplificationSettings> simulation simplification settings
    consts          <Constants> simulation's constants
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
    return
end

"""
    remove_floe_overlap!(
        floes,
        shrinking_idx,
        growing_floe_coords,
        pieces_buffer,
        max_floe_id,
        broken,
        coupling_settings,
        ridgeraft_settings,
        consts,
        rng,  
    )

Removes area/volume of overlap from floe that loses area during ridging/rafting
Inputs:
    floes               <StructArray{Floe}> list of floes
    shrinking_idx       <Int> index of floe that loses area
    growing_floe_coords <PolyVec> coordinate of floe/domain that subsumes area
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    coupling_settings   <CouplingSettings>  simulation's coupling settings
    ridgeraft_settings  <RidgeRaftSettings> simukation's ridge/raft settings
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
    ridgeraft_settings,
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
    transfer_area = (floes.area[shrinking_idx] - new_floe_area)
    transfer_vol = FT(0)
    println("in ridge")
    if transfer_area > 0.01 * floes.area[shrinking_idx]#ridgeraft_settings.min_overlap 
        println("actually ridged")
        transfer_vol = transfer_area * floes.height[shrinking_idx]
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
    end
    return transfer_vol, max_floe_id, nregions
end

"""
    floe_floe_ridge!(
        floes,
        idx1,
        idx2,
        floe2,
        overlap_area,
        ridgeraft_settings,
        simp_settings,
        consts,
        Δt,
        rng
    )
Ridge two floes, updating both in-place and returning any new floes that
resulting from the ridging event.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    idx2                <Int> index of second floe
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    ridgeraft_settings  <RidgeRaftSettings> simulation's settings for ridging
                            and rafting
    coupling_settings   <CouplingSettings> simulation's settings for coupling
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    consts              <Constants> simulation's constants
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
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
    nregions = 0
    # Determine which floe transfers mass and which gains mass
    gain_mass_idx = 0
    lose_mass_idx = 0 
    if(
        (f1_h && f2_h && rand() >= 1/(1 + (floes.height[idx1]/floes.height[idx2]))) ||
        (f1_h && !f2_h)
    )
        #=
        Either both floes are over min height and we randomly pick floe 1 to
        "subsume" the extra area, or only floe 1 is over min height and it
        gets the extra area
        =#
        gain_mass_idx = idx1
        lose_mass_idx = idx2
    elseif (f1_h && f2_h) ||  (!f1_h && f2_h)
        #=
        Either both floes are over max height and we randomly pick floe 2 to
        "subsumes" the extra area, or only floe 2 is over max height and it
        gets the extra area
        =#
        gain_mass_idx = idx2
        lose_mass_idx = idx1
    end
    if gain_mass_idx > 0 && lose_mass_idx > 0
        # Inital floe values
        ml, mg = floes.mass[lose_mass_idx], floes.mass[gain_mass_idx]
        Il, Ig, = floes.moment[lose_mass_idx], floes.moment[gain_mass_idx]
        xl, yl = floes.centroid[lose_mass_idx]
        xg, yg = floes.centroid[gain_mass_idx]
        # TODO: remove from debugging
        cl, cg = deepcopy(floes.coords[lose_mass_idx]), deepcopy(floes.coords[gain_mass_idx])
        # Ridge
        vol, max_floe_id, nregions = remove_floe_overlap!(
            floes,
            lose_mass_idx,
            floes.coords[gain_mass_idx],
            pieces_buffer,
            max_floe_id,
            broken,
            coupling_settings,
            ridgeraft_settings,
            consts,
            rng,  
        )
        if vol > 0 && nregions > 0
            add_floe_volume!(
                floes,
                gain_mass_idx,
                vol,
                simp_settings,
                consts,
                Δt,
            )
            # Conserve momentum
            first_slot = length(pieces_buffer) - nregions + 2
            if nregions < 2
                conserve_momentum_transfer_mass!(floes, lose_mass_idx, gain_mass_idx,
                    ml, mg, Il, Ig, xl, xg, yl, yg, Δt, cl, cg,
                )
            else  # floe broke, ghost floes
                conserve_momentum_transfer_mass!(floes, lose_mass_idx, gain_mass_idx,
                    ml, mg, Il, Ig, xl, xg, yl, yg, Δt, cl, cg, pieces_buffer, first_slot,
                )
            end
            if !broken[lose_mass_idx]
                update_ghost_timestep_vals!(floes, lose_mass_idx)
            end
            if !broken[gain_mass_idx]
                update_ghost_timestep_vals!(floes, gain_mass_idx)
            end
            @assert !isnan(floes.u[lose_mass_idx]) && !isnan(floes.u[gain_mass_idx]) &&
            !isnan(floes.v[lose_mass_idx]) && !isnan(floes.v[gain_mass_idx]) &&
            !isnan(floes.ξ[lose_mass_idx]) && !isnan(floes.ξ[gain_mass_idx]) &&
            !isnan(floes.p_dxdt[lose_mass_idx]) && !isnan(floes.p_dxdt[gain_mass_idx]) &&
            !isnan(floes.p_dydt[lose_mass_idx]) && !isnan(floes.p_dydt[gain_mass_idx]) &&
            !isnan(floes.p_dαdt[lose_mass_idx]) && !isnan(floes.p_dαdt[gain_mass_idx]) " post conserve $lose_mass_idx  $gain_mass_idx $(floes.u[lose_mass_idx]) $(floes.u[gain_mass_idx]) $(floes.v[lose_mass_idx]) $(floes.v[gain_mass_idx]) $(floes.ξ[lose_mass_idx]) $(floes.ξ[gain_mass_idx])"
        end
    end
    return max_floe_id
end

"""
    floe_domain_ridge!(
        floes,
        idx,
        domain_element,
        pieces_buffer,
        max_floe_id,
        broken,
        coupling_settings,
        consts,
        Δt,
        rng,
    )

Ridge a floe against a boundary or a topography element and return any excess
floes created by the ridging.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    domain_element      <AbstractDomainElement> boundary or topography element
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    coupling_settings   <CouplingSettings> simulation's settings for coupling
    ridgeraft_settings  <RidgeRaftSettings> simulation's settings for ridge/raft
    consts              <Constants> simulation's constants
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
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
    ridgeraft_settings,
    consts,
    Δt,
    rng,
)
    # Record previous values for momentum conservation
    mass_tmp = floes.mass[idx]
    moment_tmp = floes.moment[idx]
    x_tmp, y_tmp = floes.centroid[idx]
    # Ridge floe with domain element
    vol, max_floe_id, _ = remove_floe_overlap!(
        floes,
        idx,
        domain_element.coords,
        pieces_buffer,
        max_floe_id,
        broken,
        coupling_settings,
        ridgeraft_settings,
        consts,
        rng,  
    )
    if vol > 0
        # Update floe velocities to conserve momentum as domain element has no
        # momentum, but floe now has less mass if ridge was successful
        conserve_momentum_change_floe_shape!(
            mass_tmp,
            moment_tmp,
            x_tmp,
            y_tmp,
            Δt,
            LazyRow(floes, idx),
        )   
    end
    return max_floe_id
end

"""
    floe_floe_raft!(
        floes,
        idx1,
        idx2,
        pieces_buffer,
        max_floe_id,
        broken,
        coupling_settings,
        ridgeraft_settings,
        simp_settings,
        consts,
        Δt,
        rng,
    )

Raft two floes, updating both in-place and returning any new floes that
resulting from the rafting event.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    idx2                <Int> index of second floe
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    coupling_settings   <CouplingSettings> simulation's settings for coupling
    ridgeraft_settings  <RidgeRaftSettings> simulation's ridge/raft settings
    simp_settings       <SimplificationSettings> simulation's simplification
                            settings
    consts              <Constants> simulation's constants
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
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
    ridgeraft_settings,
    simp_settings,
    consts,
    Δt,
    rng,
) where {FT}
    vol = FT(0)
    nregions = 0
    # Based on height ratio, pick which floe subsumes shares area
    # Default is floe 2 subsumes mass from floe 1
    gain_mass_idx = idx2
    lose_mass_idx = idx1
    if rand(rng) >= 1/(1 + (floes.height[idx1]/floes.height[idx2]))
        # Floe 1 subsumes mass from floe 2
        gain_mass_idx = idx1
        lose_mass_idx = idx2
    end
    # Inital floe values
    ml, mg = floes.mass[lose_mass_idx], floes.mass[gain_mass_idx]
    Il, Ig, = floes.moment[lose_mass_idx], floes.moment[gain_mass_idx]
    xl, yl = floes.centroid[lose_mass_idx]
    xg, yg = floes.centroid[gain_mass_idx]
    # TODO remove post debugging
    cl, cg = floes.coords[lose_mass_idx], floes.coords[gain_mass_idx]
    # Raft
    vol, max_floe_id, nregions = remove_floe_overlap!(
        floes,
        lose_mass_idx,
        floes.coords[gain_mass_idx],
        pieces_buffer,
        max_floe_id,
        broken,
        coupling_settings,
        ridgeraft_settings,
        consts,
        rng,  
    )
    if vol > 0 && nregions > 0
        # Add extra area/volume to floe 2
        add_floe_volume!(
            floes,
            gain_mass_idx,
            vol,
            simp_settings,
            consts,
            Δt,
        )
        first_slot = length(pieces_buffer) - nregions + 2
        if nregions < 2
            conserve_momentum_transfer_mass!(floes, lose_mass_idx, gain_mass_idx,
                ml, mg, Il, Ig, xl, xg, yl, yg, Δt, cl, cg,
            )
        else  # floe broke, ghost floes
            conserve_momentum_transfer_mass!(floes, lose_mass_idx, gain_mass_idx,
                ml, mg, Il, Ig, xl, xg, yl, yg, Δt, cl, cg, pieces_buffer, first_slot,
            )
        end
        if !broken[lose_mass_idx]
            update_ghost_timestep_vals!(floes, lose_mass_idx)
        end
        if !broken[gain_mass_idx]
            update_ghost_timestep_vals!(floes, gain_mass_idx)
        end
    end
    return max_floe_id
end

"""
    floe_domain_raft!(
        floes,
        idx1,
        domain_element,
        pieces_buffer,
        max_floe_id,
        broken,
        coupling_settings,
        ridgeraft_settings,
        simp_settings,
        consts,
        Δt
        rng,
    )

Raft a floe against a boundary or a topography element and return any excess
floes created by the rafting. This is equivalent to ridging.
Inputs:
    floes               <StructArray{Floe}> floe list
    idx1                <Int> index of first floe
    domain_element      <AbstractDomainElement> boundary or topography element
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    broken              <Vector{Bool}> floe index is true if that floe has
                            broken in a previous ridge/raft interaction
    coupling_settings   <CouplingSettings> simulation's settings for coupling
    consts              <Constants> simulation's constants
    Δt                  <Int> simulation timestep in seconds
    rng                 <RandomNumberGenerator> simulation's random number
                            generator
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
    ridgeraft_settings,
    consts,
    Δt,
    rng,
) = floe_domain_ridge!(
    floes,
    idx1,
    domain_element,
    pieces_buffer,
    max_floe_id,
    broken,
    coupling_settings,
    ridgeraft_settings,
    consts,
    Δt,
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
    pieces_buffer       <StructArray{Floe}> list of new floe pieces caused by
                            breakage of floes
    n_init_floes        <Int> number of floes prior to adding ghosts
    domain              <Domain> simulation's domain
    max_floe_id         <Int> maximum floe ID before this ridging/rafting
    coupling_settings   <CouplingSettings> coupling settings
    ridgeraft_settings  <RidgeRaftSettings> ridge/raft settings
    simp_settings       <SimplificationSettings> simplification settings
    consts              <Consts> simulation's constants
    Δt                  <Int> length of timestep in seconds
    rng                 <RandomNumberGenerator> simulation's rng
Outputs:
    Updates floes post ridging and rafting and adds any new pieces to the pieces
    buffer to be made into new floes.
"""
function timestep_ridging_rafting!(
    floes,
    pieces_buffer,
    n_init_floes,
    domain,
    max_floe_id,
    coupling_settings,
    ridgeraft_settings::RidgeRaftSettings{FT},
    simp_settings,
    consts,
    Δt,
    rng = Xoshiro(),
) where {FT <: AbstractFloat}
    broken = fill(false, length(floes))
    interactions_list = zeros(Int, size(floes.interactions[1], 1))
    for i in eachindex(floes)
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

                # floes/domain overlap (not ghost interaction copied to parent)
                valid_interaction = false
                if i < j && !broken[j]
                    valid_interaction |= potential_interaction(
                        floes.centroid[i], floes.centroid[j],
                        floes.rmax[i], floes.rmax[j],
                    )
                elseif j == -1
                    valid_interaction |=
                        abs(floes.centroid[i][2] - domain.north.val) <
                        floes.rmax[i]
                elseif j == -2
                    valid_interaction |=
                        abs(floes.centroid[i][2] - domain.south.val) <
                        floes.rmax[i]
                elseif j == -3
                    valid_interaction |=
                        abs(floes.centroid[i][1] - domain.east.val) <
                        floes.rmax[i]
                elseif j == -4
                    valid_interaction |=
                        abs(floes.centroid[i][1] - domain.west.val) <
                        floes.rmax[i]
                elseif j < 0
                    valid_interaction |= potential_interaction(
                        floes.centroid[i],
                        domain.topography.centroid[-(j + 4)],
                        floes.rmax[i],
                        domain.topography.rmax[-(j + 4)],
                    )
                end
                if (valid_interaction && !(j in interactions_list) &&
                    1e-6 < floes.interactions[i][row, overlap]/min_area < 0.95
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
                            pieces_buffer,
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
                            pieces_buffer,
                            max_floe_id,
                            broken,
                            coupling_settings,
                            ridgeraft_settings,
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
                            pieces_buffer,
                            max_floe_id,
                            broken,
                            coupling_settings,
                            ridgeraft_settings,
                            consts,
                            Δt,
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
                            pieces_buffer,
                            max_floe_id,
                            broken,
                            coupling_settings,
                            ridgeraft_settings,
                            consts,
                            Δt,
                            rng,
                        )
                    end
                end
            end
        end
    end
    return max_floe_id
end
