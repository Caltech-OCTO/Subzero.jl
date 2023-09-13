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
    # TODO: Make sure to update ghosts / parents
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
    shrinking_floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    shrinking_floeidx,
    growing_floe_coords,
    vol,
    pieces_buffer,
    coupling_settings,
    consts,
    rng,  
) where {FT <: AbstractFloat}
    # TODO: Make sure to update ghosts / parents
    # Find new floe shape
    new_floe_poly = LG.difference(
        LG.Polygon(shrinking_floe.coords),
        LG.Polygon(growing_floe_coords),
    )
    new_floe_area = LG.area(new_floe_poly)
    transfer_vol = (shrinking_floe.area - new_floe_area) * shrinking_floe.height
    regions = get_polygons(new_floe_poly)
    nregions = length(regions)
    # Update floe and create new floes from ridging
    for i in 1:nregions
        new_coords = find_poly_coords(regions[i])::PolyVec{FT}
        new_poly = LG.Polygon(rmholes(new_coords))
        mass_diff =  (shrinking_floe.mass - transfer_vol * consts.ρi)
        new_mass = (LG.area(new_poly) / new_floe_area) *
            (shrinking_floe.mass - transfer_vol * consts.ρi)
        mass_tmp = shrinking_floe.mass
        moment_tmp = shrinking_floe.moment
        x_tmp, y_tmp = shrinking_floe.centroid
        if i == 1  # Update existing floe
            # Update floe based on new shape
            replace_floe!(
                shrinking_floe,
                new_poly,
                new_mass,
                coupling_settings,
                consts,
                rng,
            )
            # TODO: Update velocities and accelerations to conserve momentum
        else  # Create new floes if floe was split into multiple pieces
            push!(pieces_buffer, deepcopy_floe(shrinking_floe))
            # Update floe based on new shape
            replace_floe!(
                LazyRow(pieces_buffer, length(pieces_buffer)),
                new_poly,
                new_mass,
                coupling_settings,
                consts,
                rng,
            )
        # TODO: Update velocities and accelerations to conserve momentum
        end
    end
    return transfer_vol, nregions > 1  # broken if true
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
    broke = false
    if(
        (f1_h && f2_h && rand() >= 1/(1 + (floes.height[idx1]/floes.height[idx2]))) ||
        (f1_h && !f2_h)
    )
        #=
        Either both floes are over min height and we randomly pick floe 1 to
        "subsume" the extra area, or only floe 1 is over min height and it
        gets the extra area
        =#
        vol, broke = remove_floe_overlap!(
            LazyRow(floes, idx2),
            idx2,
            floes.coords[idx1],
            vol,
            pieces_buffer,
            coupling_settings,
            consts,
            rng,  
        )
        broken[idx2] = broke
        add_floe_volume!(
            LazyRow(floes, idx1),
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
        vol, broke = remove_floe_overlap!(
            LazyRow(floes, idx1),
            idx1,
            floes.coords[idx2],
            vol,
            pieces_buffer,
            coupling_settings,
            consts,
            rng,  
        )
        broken[idx2] = broke
        add_floe_volume!(
            LazyRow(floes, idx2),
            vol,
            simp_settings,
            consts,
            Δt,
        )
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
    floe::Union{Floe, LazyRow{Floe}},
    idx,
    domain_element::Union{
        CollisionBoundary,
        CompressionBoundary,
        TopographyElement,
    },
    pieces_buffer,
    broken,
    coupling_settings,
    simp_settings,
    consts,
    rng,
)
    _, broken[idx] = remove_floe_overlap!(
        floe,
        idx,
        domain_element.coords,
        vol,
        pieces_buffer,
        coupling_settings,
        consts,
        rng,  
    )
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
    floes::StructArray{<:Floe{FT}},
    idx1,
    idx2,
    pieces_buffer,
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
        vol, broken[idx2] = remove_floe_overlap!(
            LazyRow(floes, idx2),
            idx2,
            floes.coords[idx1],
            vol,
            pieces_buffer,
            coupling_settings,
            consts,
            rng,  
        )
        # Add extra area/volume to floe 1
        add_floe_volume!(
            LazyRow(floes, idx1),
            vol,
            simp_settings,
            consts,
            Δt,
        )
    else  # Floe 2 subsumes
        # Change shape of floe 1
        vol, broken[idx1] = remove_floe_overlap!(
            LazyRow(floes, idx1),
            idx1,
            floes.coords[idx2],
            vol,
            pieces_buffer,
            coupling_settings,
            consts,
            rng,  
        )
        # Add extra area/volume to floe 2
        add_floe_volume!(
            LazyRow(floes, idx2),
            vol,
            simp_settings,
            consts,
            Δt,
        )
    end
    return nothing
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
    floe1::Union{Floe, LazyRow{Floe}},
    idx1,
    domain_element::AbstractDomainElement,
    pieces_buffer,
    broken,
    coupling_settings,
    simp_settings,
    consts,
    rng,
) = floe_domain_ridge!(
    floe1,
    idx1,
    domain_element,
    pieces_buffer,
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
    n_init_floes,
    domain,
    ridgeraft_settings::RidgeRaftSettings{FT},
    coupling_settings,
    simp_settings,
    consts,
    Δt,
    rng = Xoshiro(),
) where {FT <: AbstractFloat}
    broken = fill(false, n_init_floes)
    pieces_list = Floe{FT}[]
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
                    (j > i || j < 0) &&  # else checked in j's interactions
                    1e-6 < floes.interactions[i][row, overlap]/min_area < 0.95 &&
                    !broken[j] && !(j in interactions_list)
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
                        floe_floe_ridge!(
                            floes,
                            i,
                            j,
                            pieces_list,
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
                        floe_floe_raft!(
                            floes,
                            i,
                            j,
                            pieces_list,
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
                        floe_domain_ridge!(
                            LazyRow(floes, i),
                            i,
                            get_domain_element(domain, j),
                            pieces_list,
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
                        floe_domain_raft!(
                            LazyRow(floes, i),
                            i,
                            get_domain_element(domain, j),
                            pieces_list,
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
    return
end
