"""
Functions needed for collisions between floes, boundaries, and topography
"""

"""
    calc_normal_force(c1, c2, region, area, ipoints, force_factor, t::Type{T} = Float64)

Calculate normal force for collision between polygons with coordinates c1 and c2 given an overlapping region,
the area of that region, their intersection points in the region, and a force factor. 
Inputs:
        c1           <PolyVec{Float64}> first polygon coordinates
        c2           <PolyVec{Float64}> second polygon coordinates
        region       <PolyVec{Float64}> coordiantes for one region of intersection between the polygons
        area         <Float> area of region
        ipoints      <Array{Float, N, 2}> Points of intersection between polygon 1 and 2 where the first column
                                          is the x-coordinates and the second column is the y-coordinates
        force_factor <Float> Spring constant equivalent for collisions
        t            <Type> Float type model is running on (Float64 or Float32)
Outputs:
        <Float> normal force of collision
        Δl <Float> mean length of distance between intersection points
"""
function calc_normal_force(c1, c2, region, area, ipoints, force_factor, ::Type{T} = Float64) where T
    force_dir = zeros(T, 2)
    coords = LG.GeoInterface.coordinates(region)::PolyVec{Float64}
    n_ipoints = size(ipoints, 1)
    # Identify which region coordinates are the intersection points (ipoints)
    verts = zeros(Int64, n_ipoints)
    dists = zeros(n_ipoints)
    for i in 1:n_ipoints
        p_i = repeat(ipoints[i, :], length(coords))
        dists[i], verts[i] = findmin([sum((c .- p_i).^2) for c in coords[1]])
    end
    dists = sqrt.(dists)
    p = coords[1][verts[findall(d -> d<1, dists)]] # maybe do rmax/1000
    m = length(p)

    # Calculate force direction
    Δl = zeros(T, 1)
    if m == 2  # Two intersection points
        Δx = p[2][1] - p[1][1]
        Δy = p[2][2] - p[1][2]
        Δl = sqrt(Δx^2 + Δy^2)
        if Δl > 0.1  # should match our scale
            force_dir = [-Δy/Δl; Δx/Δl]
        end
    elseif m != 0  # Unusual number of intersection points
        x, y = seperate_xy(coords)
        Δx = diff(x)
        xmid = (x[2:end] .+ x[1:end-1]) ./ 2
        Δy = diff(y)
        ymid = (y[2:end] .+ y[1:end-1]) ./ 2
        mag = sqrt.(Δx.^2 .+ Δy.^2)  # vector magniture
        uvec = [-Δy./mag Δx./mag]  # unit vector
        xt = xmid.+uvec[:, 1]./100
        yt = ymid+uvec[:, 2]./100  # should match our scale
        in_idx = inpoly2(hcat(xt, yt), hcat(x, y))
        uvec[in_idx[:, 1] .|  in_idx[:, 2], :] *= -1
        Fn = -force_factor * (mag * ones(1, 2)) .* uvec
        dmin_lst = calc_point_poly_dist(xmid, ymid, c1)
        on_idx = findall(d->abs(d)<1e-8, dmin_lst)
        if 0 < length(on_idx) < length(dmin_lst)
            Δl = mean(mag[on_idx])
            if Δl > 0.1
                Fn_tot = sum(Fn[on_idx, :], dims = 1)
                force_dir = Fn_tot'/sqrt(Fn_tot*Fn_tot')
            end
        end
    end
    # Check if direction of the force desceases overlap, else negate direction
    if Δl > 0.1
        c1new = [[c .+ vec(force_dir) for c in c1[1]]]
        # Floe/boudary intersection after being moved in force direction
        new_inter_floe = LG.intersection(LG.Polygon(c1new), LG.Polygon(c2))
        # See if the area of overlap has increased in corresponding region
        for new_region in LG.getGeometries(new_inter_floe)
            if LG.intersects(new_region, region) && LG.area(new_region)/area > 1
                force_dir *= -1 
            end
        end
    end
    return (force_dir * area * force_factor)', Δl
end

"""
    calc_elastic_forces(floe1, floe2, regions, region_areas, force_factor, consts, t::Type{T} = Float64)

Calculate normal forces, the point the force is applied, and the overlap area of regions created from floe collisions 
Inputs:
        floe1           <Floe> first floe in collision
        floe2           <Floe> second floe in collision
        regions         <Vector{LibGEOS.Polygon}> polygon regions of overlap during collision
        region_areas    <Vector{Float}> area of each polygon in regions
        force_factor    <Float> Spring constant equivalent for collisions
        t               <Type> Float type model is running on (Float64 or Float32)
Outputs:
        force   <Array{Float, n, 2}> normal forces on each of the n regions greater than a minimum area
        fpoint  <Array{Float, n, 2}> point force is applied on each of the n regions greater than a minimum area
        overlap <Array{Float, n, 2}> area of each of the n regions greater than a minimum area
        Δl      <Float> mean length of distance between intersection points
"""
function calc_elastic_forces(c1, c2, regions, region_areas, force_factor, ::Type{T} = Float64) where T
    ipoints = intersect_lines(c1, c2)  # Intersection points
    if isempty(ipoints) || size(ipoints,2) < 2  # No overlap points
        return zeros(T, 1, 2), zeros(T, 1, 2), zeros(T, 1)  # Force, contact points, overlap area
    else
        # Find overlapping regions greater than minumum area
        n1 = length(c1[1]) - 1
        n2 = length(c2[1]) - 1
        min_area = min(n1, n2) * 100 / 1.75
        regions = regions[region_areas .> min_area]
        region_areas = region_areas[region_areas .> min_area]
        overlap = region_areas
        ncontact = length(regions)
        # Calculate forces for each remaining region
        force = zeros(T, ncontact, 2)
        fpoint = zeros(T, ncontact, 2)
        Δl_lst = zeros(T, ncontact)
        for k in 1:ncontact
            normal_force = zeros(T, 1, 2)
            if region_areas[k] != 0
                cx, cy = LG.GeoInterface.coordinates(LG.centroid(regions[k]))::Vector{Float64}
                fpoint[k, :] = [cx, cy]
                normal_force, Δl = calc_normal_force(c1, c2, regions[k], region_areas[k], ipoints, force_factor, T)
            end
            force[k, :] = normal_force
            Δl_lst[k] = Δl
        end
        return force, fpoint, overlap, Δl_lst
    end
end

"""
    calc_friction_forces(v1, v2, fpoint, nforce, Δl, consts)

Calculate frictional force for collision between two floes.
Input:
        v1   <Array{Float, N, 2}> Matrix of floe speeds for first floe in collision at each point
                                  a force is applied to the floe - row is [u v] and one row per each N point
        v2   <Array{Float, N, 2}> vector of floe speeds for second floe or boundary in collision - see v1 for form
        fpoint  <Array{Float, N, 2}> x,y-coordinates of the point the force is applied on floe overlap region
        normal  <Array{Float, N, 2}> x,y normal force applied on fpoint on floe overlap region
        Δl      <Float> mean length of distance between intersection points
        consts  <Constants> model constants needed for calculations
Outputs:
        <Float> frictional/tangential force of the collision
"""
function calc_friction_forces(v1, v2, normal, Δl, consts, Δt, ::Type{T} = Float64) where T
    force = zeros(T, size(v1, 1), 2)
    G = consts.E/(2*(1+consts.ν))  # Sheer modulus
    # Difference in velocities between floes in x and y direction
    vdiff = v1 .- v2
    # Friction forces for each vector
    for i in axes(vdiff, 1)
        v = vdiff[i, :]
        n = normal[i, :]
        vnorm = norm(v)
        force_dir = maximum(abs.(v)) == 0 ? zeros(T, 2) : v/vnorm
        friction = G * Δl[i] * Δt * vnorm * -dot(force_dir, v) * force_dir
        if norm(friction) > consts.μ*norm(n)
            friction = -consts.μ*norm(n)*force_dir
        end
        force[i, :] = friction
    end
    return force
end

"""
    floe_floe_interaction!(ifloe, i, jfloe, j, nfloes, consts, Δt, t::Type{T} = Float64)

If the two floes interact, update floe i's interactions accordingly. Floe j is not update here so that the function can be parallelized.
Inputs:
        ifloe   <Floe> first floe in potential interaction
        i       <Int> index of ifloe in model's list of floes 
        jfloe   <Floe> second floe in potential interaction
        j       <Int> index of jfloe in model's list of floes 
        nfloes  <Int> number of non-ghost floes in the simulation this timestep
        consts  <Constants> model constants needed for calculations
        Δt      <Int> Simulation's current timestep
        t       <Type> Float type model is running on (Float64 or Float32)
Outputs:
        remove   <Int> index of floe to remove from simulation due to overlap (will be transfered to over floe in collision)
                        - if 0 then no floes to be removed from this collision
        transfer <Int> index of floe to transfer mass to if the other floe is to be removed.
                        - if 0 then no floes to be transfered from this collision
        If ifloe interacts with jfloe, only ifloe's interactions field is updated with the details of each region of overlap.
        The interactions field will have the following form for each region of overlap with the boundary:
            [Inf, xforce, yforce, xfpoints, yfpoints, overlaps] where the xforce and yforce are the forces,
            xfpoints and yfpoints are the location of the force and overlaps is the overlap between the floe and boundary.
        The overlaps field is also added to the floe's overarea field that describes the total overlapping area at any timestep. 
"""
function floe_floe_interaction!(ifloe, i, jfloe, j, nfloes, consts, Δt, ::Type{T} = Float64) where T
    remove = Int(0)
    transfer = Int(0)
    ifloe_poly = LG.Polygon(ifloe.coords)
    jfloe_poly = LG.Polygon(jfloe.coords)
    if LG.intersects(ifloe_poly, jfloe_poly)  # Check if floes intersect
        inter_floe = LG.intersection(ifloe_poly, jfloe_poly)
        inter_regions = LG.getGeometries(inter_floe)
        region_areas = [LG.area(poly) for poly in inter_regions]::Vector{Float64}
        total_area = sum(region_areas)
        # Floes overlap too much - remove floe or transfer floe mass to other floe
        if total_area/ifloe.area > 0.55
            if i <= nfloes  # If i is not a ghost floe
                remove = i
                transfer = j  # Will transfer mass to jfloe
            elseif j <= nfloes  # If j is not a ghost floe
                remove = j  # Will transfer mass to ifloe - can't be updated here due to parallelization
            end
        elseif total_area/jfloe.area > 0.55
            if j <= nfloes  # If j is not a ghost floe
                remove = j  # Will transfer mass to ifloe - can't be updated here due to parallelization
            end
        else
            # Constant needed to calculate force
            ih = ifloe.height
            ir = sqrt(ifloe.area)
            jh = jfloe.height
            jr = sqrt(jfloe.area)
            force_factor = if ir>1e5 || jr>1e5
                consts.E*min(ih, jh)/min(ir, jr)
            else
                consts.E*(ih*jh)/(ih*jr+jh*ir)
            end
            # Calculate normal forces, force points, and overlap areas
            normal_forces, fpoints, overlaps, Δl = calc_elastic_forces(ifloe.coords, jfloe.coords,
                                        inter_regions, region_areas, force_factor, T)
            # Calculate frictional forces at each force point - based on velocities at force points
            np = size(fpoints, 1)
            if np > 0
                iv = repeat([ifloe.u ifloe.v], outer = np) .+ ifloe.ξ*(fpoints .- repeat(ifloe.centroid', outer = np)) 
                jv = repeat([jfloe.u jfloe.v], outer = np) .+ jfloe.ξ*(fpoints .- repeat(jfloe.centroid', outer = np))
                friction_forces = calc_friction_forces(iv, jv, normal_forces, Δl, consts, Δt, T)
                # Calculate total forces and update ifloe's interactions
                forces = normal_forces .+ friction_forces
                if sum(abs.(forces)) != 0
                    ifloe.interactions = [ifloe.interactions; fill(j, np) forces fpoints zeros(np) overlaps]
                    ifloe.overarea += sum(overlaps)
                end
            end
        end
    end
    return remove, transfer
end

"""
    floe_domain_element_interaction!(floe, boundary, _, _, _)

If given floe insersects with an open boundary, the floe is set to be removed from the simulation.
Inputs:
        floe            <Floe> floe interacting with boundary
        boundary        <OpenBoundary> coordinates of boundary
        _               <Constants> model constants needed in other methods of this function - not needed here
        _               <Int> current simulation timestep - not needed here
        _               <Type> Float type model is running on (Float64 or Float32) - not needed here
Output:
        None. If floe is interacting with the boundary, floe.alive field is set to 0. Else, nothing is changed. 
"""
function floe_domain_element_interaction!(floe, boundary::OpenBoundary, consts, Δt, ::Type{T} = Float64) where T
    floe_poly = LG.Polygon(floe.coords)
    bounds_poly = LG.Polygon(boundary.coords)
    # Check if the floe and boundary actually overlap
    if LG.intersects(floe_poly, bounds_poly)
        floe.alive = 0
    end
    return
end

"""
    floe_domain_element_interaction!(floe, ::PeriodicBoundary, consts, Δt, ::Type{T} = Float64)

If a given floe intersects with a periodic boundary, nothing happens at this point. Periodic floes pass through boundaries
using ghost floes.
Inputs:
        None are used. 
Output:
        None. This function does not do anyting. 
"""
function floe_domain_element_interaction!(floe, ::PeriodicBoundary, consts, Δt, ::Type{T} = Float64) where T
    return
end

"""
    normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{North, <:AbstractFloat}, ::Type{T} = Float64)

Zero-out forces that point in direction not perpendicular to North or South boundary wall.
Inputs:
        force       <Array{Float, n, 2}> normal forces on each of the n regions greater than a minimum area
        fpoint      <Array{Float, n, 2}> point force is applied on each of the n regions greater than a minimum area
        boundary    <AbstractBoundary{North, <:AbstractFloat}> domain's northern boundary 
                    <Type> Float type model is running on (Float64 or Float32) - not needed here
Outputs: None. All forces in the x direction set to 0 if the point the force is applied is the northern or southern boundary value.
"""
function normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{North, <:AbstractFloat}, ::Type{T} = Float64) where T
    forces[fpoints[:, 2] .>= boundary.val, 1] .= T(0.0)
    return
end

function normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{South, <:AbstractFloat}, ::Type{T} = Float64) where T
        forces[fpoints[:, 2] .<= boundary.val, 1] .= T(0.0)
        return
    end

"""
    normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{Union{East, West}, <:AbstractFloat}, ::Type{T} = Float64)

Zero-out forces that point in direction not perpendicular to East or West boundary wall.
Inputs:
        force       <Array{Float, n, 2}> normal forces on each of the n regions greater than a minimum area
        fpoint      <Array{Float, n, 2}> point force is applied on each of the n regions greater than a minimum area
        boundary    <AbstractBoundary{East, <:AbstractFloat}> domain's southern boundary 
                    <Type> Float type model is running on (Float64 or Float32) - not needed here
Outputs: None. All forces in the y direction set to 0 if the point the force is applied is the eastern or western boundary value.
"""
function normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{East, <:AbstractFloat}, ::Type{T} = Float64) where T
    forces[fpoints[:, 1] .>= boundary.val, 2] .= T(0.0)
    return
end

function normal_direction_correct!(forces, fpoints, boundary::Union{AbstractBoundary{East, <:AbstractFloat},
    AbstractBoundary{West, <:AbstractFloat}}, ::Type{T} = Float64) where T
forces[fpoints[:, 1] .<= boundary.val, 2] .= T(0.0)
return
end

"""
    normal_direction_correct!(forces, fpoints, ::TopographyElement, ::Type{T} = Float64)

No forces should be zero-ed out in collidions with topography elements. 
Inputs:
        None used.
Outputs:
        None.
"""
function normal_direction_correct!(forces, fpoints, ::TopographyElement, ::Type{T} = Float64) where T
    return
end

"""
    floe_domain_element_interaction!(floe, element, consts, Δt, t::Type{T} = Float64)

If floe intersects with given element (either collision boundary or topography element), floe interactions field and overarea field are updated.
Inputs:
        floe            <Floe> floe interacting with element
        element         <Union{CollisionBoundary, TopographyElement}> coordinates of element
        consts          <Constants> model constants needed for calculations
        Δt              <Int> current simulation timestep
                        <Type> Float type model is running on (Float64 or Float32)
Outputs:
        None. If floe interacts, the floe's interactions field is updated with the details of each region of overlap.
        The interactions field will have the following form for each region of overlap with the element:
            [Inf, xforce, yforce, xfpoints, yfpoints, overlaps] where the xforce and yforce are the forces,
            xfpoints and yfpoints are the location of the force and overlaps is the overlap between the floe and element.
        The overlaps field is also added to the floe's overarea field that describes the total overlapping area at any timestep.
"""
function floe_domain_element_interaction!(floe, element::Union{CollisionBoundary, TopographyElement}, consts, Δt, ::Type{T} = Float64) where T
    floe_poly = LG.Polygon(floe.coords)
    bounds_poly = LG.Polygon(element.coords)
    # Check if the floe and element actually overlap
    if LG.intersects(floe_poly, bounds_poly)
        inter_floe = LG.intersection(floe_poly, bounds_poly)
        inter_regions = LG.getGeometries(inter_floe)
        region_areas = [LG.area(poly) for poly in inter_regions]::Vector{Float64}
        # Regions overlap too much
        if maximum(region_areas)/floe.area > 0.75
            floe.alive = 0
        else
            # Constant needed for force calculations
            force_factor = consts.E * floe.height / sqrt(floe.area)
            # Calculate normal forces, force points, and overlap areas
            normal_forces, fpoints, overlaps, Δl =  calc_elastic_forces(floe.coords, element.coords,
                                            inter_regions, region_areas, force_factor, T)
            normal_direction_correct!(normal_forces, fpoints, element, T)
            # Calculate frictional forces at each force point - based on velocities at force points
            np = size(fpoints, 1)
            vfloe = repeat([floe.u floe.v], outer = np) .+ floe.ξ*(fpoints .- repeat(floe.centroid', outer = np)) 
            vbound = repeat(zeros(T, 1, 2), outer = np)
            friction_forces = calc_friction_forces(vfloe, vbound, normal_forces, Δl, consts, Δt, T)
            # Calculate total forces and update ifloe's interactions
            forces = normal_forces .+ friction_forces
            if sum(abs.(forces)) != 0
                floe.interactions = [floe.interactions; fill(Inf, np) forces fpoints zeros(np) overlaps]
                floe.overarea += sum(overlaps)
            end
        end
    end
    return
end

"""
    floe_domain_interaction!(floe, domain::DT, consts, t::Type{T} = Float64)

If the floe interacts with the domain, update the floe accordingly. Dispatches on different boundary types within the domain.
Inputs:
        floe        <Floe> floe interacting with boundary
        domain      <Domain> model domain
        consts      <Constants> model constants needed for calculations
        Δt          <Int> current simulation timestep
        t           <Type> Float type model is running on (Float64 or Float32)
Outputs:
        None. Floe is updated according to which boundaries it interacts with and the types of those boundaries. 
"""
function floe_domain_interaction!(floe, domain::Domain, consts, Δt, ::Type{T} = Float64) where {T}
    centroid = floe.centroid
    rmax = floe.rmax
    nbound = domain.north
    sbound = domain.south
    ebound = domain.east
    wbound = domain.west

    if centroid[2] + rmax > nbound.val
        floe_domain_element_interaction!(floe, nbound, consts, Δt)
    end
    if centroid[2] - rmax < sbound.val
        floe_domain_element_interaction!(floe, sbound, consts, Δt)
    end
    if centroid[1] + rmax > ebound.val
        floe_domain_element_interaction!(floe, ebound, consts, Δt)
    end
    if centroid[1] - rmax < wbound.val
        floe_domain_element_interaction!(floe, wbound, consts, Δt)
    end

    for topo_element in domain.topography
        if sum((topo_element.centroid .- floe.centroid).^2) < (topo_element.rmax + floe.rmax)^2
            floe_domain_element_interaction!(floe, topo_element, consts, Δt)
        end
    end

    return
end

"""
    calc_torque!(floe, t::Type{T} = Float64)

Calculate a floe's torque based on the interactions.
Inputs:
        floe  <Floe> floe in model
              <Type> Float type model is running on (Float64 or Float32)
Outputs:
        None. Floe's interactions field is updated with calculated torque.
"""
function calc_torque!(floe, ::Type{T} = Float64) where T
    inters = floe.interactions
    if !isempty(inters)
        dir = [inters[:, xpoint] .- floe.centroid[1] inters[:, ypoint] .- floe.centroid[2] zeros(T, size(inters, 1))]
        frc = [inters[:, xforce] inters[:, yforce] zeros(T, size(inters, 1))]
        for i in axes(dir, 1)
            idir = vec(dir[i, :])
            ifrc = vec(frc[i, :])
            itorque = cross(idir, ifrc)
            floe.interactions[i, torque] = itorque[3]
        end
    end
end

function timestep_collisions!(floes, n_init_floes, domain, remove, transfer, consts, Δt, ::Type{T} = Float64) where T
    collide_pairs = Set{Tuple{Int, Int}}()
    # floe-floe collisions for floes i and j where i<j
    for i in eachindex(floes)
        ifloe = floes[i]
        iid = abs(ifloe.id)
        # reset collision values
        ifloe.collision_force = zeros(T, 1, 2)
        ifloe.collision_trq = T(0.0)
        ifloe.interactions = zeros(T, 0, 7)
        for j in i+1:length(floes)
            jid = abs(floes[j].id)
            # Confirm that this collision hasn't already occured with parent/ghost floes using ids
            # If it hasn't occured, check if floes are in close proximity
            if (jid != iid) && !((jid, iid) in collide_pairs) && !((iid, jid) in collide_pairs) &&
            (sum((ifloe.centroid .- floes[j].centroid).^2) < (ifloe.rmax + floes[j].rmax)^2)
                push!(collide_pairs, (iid, jid))  # add ids to collide_pairs list so collision isn't repeated with ghosts
                iremove, itransfer = floe_floe_interaction!(ifloe, i, floes[j], j, n_init_floes, consts, Δt)
                if iremove != 0 || itransfer != 0
                    remove[i] = iremove
                    transfer[i] = itransfer
                end
            end
        end
        floe_domain_interaction!(ifloe, domain, consts, Δt)
        floes[i] = ifloe
    end
    # Update floes not directly calculated above where i>j - can't be parallelized
    for i in eachindex(floes)
        # if we are removing floe j = remove[i] floe j's mass must be transfered to floe i
        if i <= n_init_floes && remove[i] > 0 && remove[i] != i
            transfer[remove[i]] = i
        end
        ij_inters = floes[i].interactions
        if !isempty(ij_inters)
            for inter_idx in axes(ij_inters, 1)  # Loop over each interaction with Floe i
                j = ij_inters[inter_idx, floeidx]  # Index of floe to update in model floe list
                if j <= length(floes) && j > i
                    jidx = Int(j)
                    floes.interactions[jidx] = [floes.interactions[jidx]; i -ij_inters[inter_idx, xforce] -ij_inters[inter_idx, yforce] #=
                                             =# ij_inters[inter_idx, xpoint] ij_inters[inter_idx, ypoint] T(0.0) ij_inters[inter_idx, overlap]]
                    floes.overarea[jidx] += ij_inters[inter_idx, overlap]
                end
            end
        end
    end
    # reverse loop so ghost torques and foreces are calculated before parents
    for i in range(1, n_init_floes)
        ifloe = floes[i]
        for g in ifloe.ghosts
            g_inters = floes.interactions[g]
            g_inters[:, xpoint] .-= (floes.centroid[g][1] - ifloe.centroid[1])
            g_inters[:, ypoint] .-= (floes.centroid[g][2] - ifloe.centroid[2])
            ifloe.interactions = [ifloe.interactions; g_inters]
        end
        calc_torque!(ifloe)
        ifloe.collision_force[1] += sum(ifloe.interactions[:, xforce])
        ifloe.collision_force[2] += sum(ifloe.interactions[:, yforce])
        ifloe.collision_trq += sum(ifloe.interactions[:, torque])
        floes[i] = ifloe
    end
    return remove, transfer
end

"""
    ghosts_on_bounds(element, ghosts, boundary, trans_vec)

If the given element intersects with the boundary, add ghosts of the element and 
any of its existing ghosts. 
Inputs:
        element     <StructArray{Floe} or StructArray{TopographyElement}> given element
        ghosts      <StructArray{Floe} or StructArray{TopographyElement}> current ghosts of element
        boundary    <PeriodicBoundary> boundary to translate element through
        trans_vec   <Matrix{Float}> 1x2 matrix of form [x y] to translate element through the boundary
Outputs:
        New ghosts created by the given element, or its current ghosts, passing through the given boundary. 
"""
function ghosts_on_bounds(element, ghosts, boundary, trans_vec)
    new_ghosts = StructArray(Vector{typeof(element)}())
    if LG.intersects(LG.Polygon(element.coords), LG.Polygon(boundary.coords))
        # ghosts of existing ghosts and original element
        append!(new_ghosts, deepcopy.(ghosts))
        push!(new_ghosts, deepcopy(element))
        for i in eachindex(new_ghosts)
            new_ghosts.coords[i] = translate(new_ghosts.coords[i], trans_vec)
            new_ghosts.centroid[i] .+= trans_vec
        end
    end
    return new_ghosts
end

"""
    find_ghosts(elem, current_ghosts, ebound::PeriodicBoundary{East, <:AbstractFloat},
                wbound::PeriodicBoundary{West, <:AbstractFloat}, ::Type{T} = Float64)

Find ghosts of given element and its known ghosts through an eastern or western periodic boundary.
If element's centroid isn't within the domain in the east/west direction, swap it with its ghost since
the ghost's centroid must then be within the domain. 
Inputs:
        elem            <StructArray{Floe} or StructArray{TopographyElement}> given element
        current_ghosts  <StructArray{Floe} or StructArray{TopographyElement}> current ghosts of element
        eboundary        <PeriodicBoundary{East, Float}> domain's eastern boundary
        wboundary        <PeriodicBoundary{West, Float}> domain's western boundary
Outputs:
        Return "primary" element, which has its centroid within the domain in the east/west direction,
        and all of its ghosts in the east/west direction, including ghosts of previously existing ghosts.
"""
function find_ghosts(elem, current_ghosts, ebound::PeriodicBoundary{East, <:AbstractFloat}, wbound::PeriodicBoundary{West, <:AbstractFloat},
::Type{T} = Float64) where T
    Lx = ebound.val - wbound.val
    new_ghosts =
        if elem.centroid[1] - elem.rmax < wbound.val  # passing through western boundary
            ghosts_on_bounds(elem, current_ghosts, wbound, [Lx, T(0)])
        elseif (elem.centroid[1] + elem.rmax > ebound.val) # passing through eastern boundary
            ghosts_on_bounds(elem, current_ghosts, ebound, [-Lx, T(0)])
        else
            StructArray(Vector{typeof(elem)}())
        end
    # if element centroid isn't in domain in east/west direction, swap with its ghost
    if !isempty(new_ghosts) && ((elem.centroid[1] < wbound.val) || (ebound.val < elem.centroid[1]))
        elem, new_ghosts[end] = new_ghosts[end], elem
    end
    return elem, new_ghosts
end

"""
    find_ghosts(elem, current_ghosts, nbound::PeriodicBoundary{North, <:AbstractFloat},
                sbound::PeriodicBoundary{South, <:AbstractFloat}, ::Type{T} = Float64)

Find ghosts of given element and its known ghosts through an northern or southern periodic boundary.
If element's centroid isn't within the domain in the north/south direction, swap it with its ghost since
the ghost's centroid must then be within the domain. 
Inputs:
        elem            <StructArray{Floe} or StructArray{TopographyElement}> given element
        current_ghosts  <StructArray{Floe} or StructArray{TopographyElement}> current ghosts of element
        nboundary        <PeriodicBoundary{North, Float}> domain's northern boundary
        sboundary        <PeriodicBoundary{South, Float}> domain's southern boundary
Outputs:
        Return "primary" element, which has its centroid within the domain in the north/south direction,
        and all of its ghosts in the north/south direction, including ghosts of previously existing ghosts.
"""
function find_ghosts(elem, current_ghosts, nbound::PeriodicBoundary{North, <:AbstractFloat}, sbound::PeriodicBoundary{South, <:AbstractFloat},
::Type{T} = Float64) where T
    Ly =  nbound.val - sbound.val
    new_ghosts = 
        if (elem.centroid[2] - elem.rmax < sbound.val)  # passing through southern boundary
            ghosts_on_bounds(elem, current_ghosts, sbound, [T(0), Ly])
        elseif (elem.centroid[2] + elem.rmax > nbound.val)  # passing through northern boundary
            ghosts_on_bounds(elem, current_ghosts, nbound, [T(0), -Ly])
        else
            StructArray(Vector{typeof(elem)}())
        end
    # if element centroid isn't in domain in north/south direction, swap with its ghost
    if !isempty(new_ghosts) && ((elem.centroid[2] < sbound.val) || (nbound.val < elem.centroid[2]))
        elem, new_ghosts[end] = new_ghosts[end], elem
    end
    return elem, new_ghosts
end

"""
    add_floe_ghosts!(floes::StructArray{Floe{T}}, max_boundary, min_boundary) where {T <: AbstractFloat}

Add ghosts of all of the given floes passing through the two given boundaries to the list of floes.
Inputs:
        floes           <StructArray{Floe{T}}> list of floes to find ghosts for
        max_boundary    <PeriodicBoundary> either northern or eastern boundary  of domain
        min_boundary    <PeriodicBoundary> either southern or western boundary of domain
Outputs:
        None. Ghosts of floes are added to floe list. 
"""
function add_floe_ghosts!(floes::StructArray{Floe{T}}, max_boundary, min_boundary) where {T <: AbstractFloat}
    nfloes = length(floes)
    for i in eachindex(floes)  # uses initial length of floes so we can append to list
        f = floes[i]
        if f.alive == 1 && f.id >= 0
            f, new_ghosts = find_ghosts(f, floes[f.ghosts], max_boundary, min_boundary, T)
            if !isempty(new_ghosts)
                new_ghosts.id .= -abs.(new_ghosts.id)  # ghosts have negative index of parent floe
                empty!.(new_ghosts.ghosts)  # remove ghost floes ghosts as these were added to parent
                append!(floes, new_ghosts)  # add ghosts to floe list
                append!(f.ghosts, nfloes+1:nfloes+length(new_ghosts))  # index of ghosts floes saved in parent
                nfloes += length(new_ghosts)
                floes[i] = f
            end
        end
    end
    return
end

"""
    add_ghosts!(elems, ::Domain{FT, <:NonPeriodicBoundary, <:NonPeriodicBoundary, <:NonPeriodicBoundary, <:NonPeriodicBoundary})

When there are no periodic boundaries, no ghosts should be added.
Inputs:
        None are used. 
Outputs:
        None. 
"""
function add_ghosts!(elems, ::Domain{FT, <:NonPeriodicBoundary, <:NonPeriodicBoundary, <:NonPeriodicBoundary, <:NonPeriodicBoundary}) where {FT<:AbstractFloat}
    return
end

"""
    add_ghosts!(elems, domain::Domain{FT, <:PeriodicBoundary, <:PeriodicBoundary, <:NonPeriodicBoundary, <:NonPeriodicBoundary})

Add ghosts for elements that pass through the northern or southern boundaries.
Inputs:
        elems   <StructArray{Floe} or StructArray{TopographyElement}> list of elements to add ghosts to
        domain  <Domain{Float, PeriodicBoundary, PeriodicBoundary,
                               NonPeriodicBoundary, NonPeriodicBoundary}> domain with northern and southern periodic boundaries
Outputs:
        None. Ghosts are added to list of elements.
"""
function add_ghosts!(elems, domain::Domain{FT, <:PeriodicBoundary, <:PeriodicBoundary, <:NonPeriodicBoundary, <:NonPeriodicBoundary}) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.north, domain.south)
    return
end

"""
    add_ghosts!(elems, domain::Domain{FT, <:NonPeriodicBoundary, <:NonPeriodicBoundary, <:PeriodicBoundary, <:PeriodicBoundary})

Add ghosts for elements that pass through the eastern or western boundaries. 
Inputs:
        elems   <StructArray{Floe} or StructArray{TopographyElement}> list of elements to add ghosts to
        domain  <Domain{Float, NonPeriodicBoundary, NonPeriodicBoundary,
                               PeriodicBoundary, PeriodicBoundary}> domain with eastern and western periodic boundaries 
Outputs:
        None. Ghosts are added to list of elements.
"""
function add_ghosts!(elems, domain::Domain{FT, <:NonPeriodicBoundary, <:NonPeriodicBoundary, <:PeriodicBoundary, <:PeriodicBoundary}) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.east, domain.west)
    return
end

"""
    add_ghosts!(elems, domain::Domain{FT, <:PeriodicBoundary, <:PeriodicBoundary, <:PeriodicBoundary, <:PeriodicBoundary})

Add ghosts for elements that pass through any of the boundaries. 
Inputs:
        elems   <StructArray{Floe} or StructArray{TopographyElement}> list of elements to add ghosts to
        domain  <Domain{Float, PeriodicBoundary, PeriodicBoundary,
                               PeriodicBoundary, PeriodicBoundary}> domain with all boundaries
Outputs:
        None. Ghosts are added to list of elements.
"""
function add_ghosts!(elems, domain::Domain{FT, <:PeriodicBoundary, <:PeriodicBoundary, <:PeriodicBoundary, <:PeriodicBoundary}) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.east, domain.west)
    add_floe_ghosts!(elems, domain.north, domain.south)
    return
end
