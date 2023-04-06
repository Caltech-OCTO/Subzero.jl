"""
Functions needed for collisions between floes, boundaries, and topography
"""

"""
    calc_normal_force(
        c1,
        c2,
        region,
        area,
        ipoints,
        force_factor,
    )

Calculate normal force for collision between polygons with coordinates c1 and c2
given an overlapping region, the area of that region, their intersection points
in the region, and a force factor. 
Inputs:
    c1           <PolyVec{Float64}> first polygon coordinates
    c2           <PolyVec{Float64}> second polygon coordinates
    region       <PolyVec{Float64}> coordiantes for one region of intersection
                    between the polygons
    area         <Float> area of region
    ipoints      <Array{Float, N, 2}> Points of intersection between polygon 1
                    and 2 where the first column is the x-coordinates and the
                    second column is the y-coordinates
    force_factor <Float> Spring constant equivalent for collisions
    t            <Type> Float type model is running on (Float64 or Float32)
Outputs:
        <Float> normal force of collision
        Δl <Float> mean length of distance between intersection points
"""
function calc_normal_force(
    c1,
    c2,
    region,
    area::FT,
    ipoints,
    force_factor::FT,
) where {FT<:AbstractFloat}
    force_dir = zeros(FT, 2)
    coords = find_poly_coords(region)
    n_ipoints = size(ipoints, 1)
    # Identify which region coordinates are the intersection points (ipoints)
    verts = zeros(Int64, n_ipoints)
    dists = zeros(FT, n_ipoints)
    for i in 1:n_ipoints
        p_i = repeat(ipoints[i, :], length(coords))
        dists[i], verts[i] = findmin([sum((c .- p_i).^2) for c in coords[1]])
    end
    dists = sqrt.(dists)
    p = coords[1][verts[findall(d -> d<1, dists)]] # maybe do rmax/1000
    m = length(p)

    # Calculate force direction
    Δl = FT(0)
    if m == 2  # Two intersection points
        Δx = p[2][1] - p[1][1]
        Δy = p[2][2] - p[1][2]
        Δl = sqrt(Δx^2 + Δy^2)
        if Δl > 0.1  # should match our scale
            force_dir = [-Δy/Δl; Δx/Δl]
        end
    elseif m != 0  # Unusual number of intersection points
        x, y = separate_xy(coords)
        Δx = diff(x)
        xmid = (x[2:end] .+ x[1:end-1]) ./ 2
        Δy = diff(y)
        ymid = (y[2:end] .+ y[1:end-1]) ./ 2
        mag = sqrt.(Δx.^2 .+ Δy.^2)  # vector magniture
        uvec = [-Δy./mag Δx./mag]  # unit vector
        xt = xmid.+uvec[:, 1]./100
        yt = ymid+uvec[:, 2]./100  # should match our scale
        in_idx = points_in_poly(hcat(xt, yt), coords)
        uvec[in_idx, :] *= -1
        Fn = -force_factor * (mag * ones(FT, 1, 2)) .* uvec
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
    calc_elastic_forces(
        floe1,
        floe2,
        regions,
        region_areas,
        force_factor,
        consts,
    )

Calculate normal forces, the point the force is applied, and the overlap area of
regions created from floe collisions 
Inputs:
    c1              <PolyVec> first floe's coordinates in collision
    c2              <PolyVec> second floe's coordinates in collision
    regions         <Vector{LibGEOS.Polygon}> polygon regions of overlap during
                        collision
    region_areas    <Vector{Float}> area of each polygon in regions
    force_factor    <Float> Spring constant equivalent for collisions
Outputs:
    force   <Array{Float, n, 2}> normal forces on each of the n regions greater
                than a minimum area
    fpoint  <Array{Float, n, 2}> point force is applied on each of the n regions
                greater than a minimum area
    overlap <Array{Float, n, 2}> area of each of the n regions greater than a
                minimum area
    Δl      <Float> mean length of distance between intersection points
"""
function calc_elastic_forces(
    c1,
    c2,
    regions,
    region_areas::Vector{FT},
    force_factor::FT,
) where {FT<:AbstractFloat}
    ipoints = intersect_lines(c1, c2)  # Intersection points
    if isempty(ipoints) || size(ipoints,2) < 2  # No overlap points
         # Force, contact points, overlap area all 0s
        return zeros(FT, 1, 2), zeros(FT, 1, 2), zeros(FT, 1)
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
        force = zeros(FT, ncontact, 2)
        fpoint = zeros(FT, ncontact, 2)
        Δl_lst = zeros(FT, ncontact)
        for k in 1:ncontact
            normal_force = zeros(FT, 1, 2)
            if region_areas[k] != 0
                cx, cy = find_poly_centroid(regions[k])::Vector{Float64}
                fpoint[k, :] = [cx, cy]
                normal_force, Δl = calc_normal_force(
                    c1,
                    c2,
                    regions[k],
                    region_areas[k],
                    ipoints,
                    force_factor
                )
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
    v1   <Array{Float, N, 2}> Matrix of floe speeds for first floe in collision
            at each point a force is applied to the floe - row is [u v] and one
            row per each N point
    v2   <Array{Float, N, 2}> vector of floe speeds for second floe or boundary
            in collision - see v1 for form
    fpoint  <Array{Float, N, 2}> x,y-coordinates of the point the force is
                applied on floe overlap region
    normal  <Array{Float, N, 2}> x,y normal force applied on fpoint on floe
                overlap region
    Δl      <Float> mean length of distance between intersection points
    consts  <Constants> model constants needed for calculations
Outputs:
        <Float> frictional/tangential force of the collision
"""
function calc_friction_forces(
    v1,
    v2,
    normal::Matrix{FT},
    Δl,
    consts,
    Δt,
) where {FT}
    force = zeros(FT, size(v1, 1), 2)
    G = consts.E/(2*(1+consts.ν))  # Sheer modulus
    # Difference in velocities between floes in x and y direction
    vdiff = v1 .- v2
    # Friction forces for each vector
    for i in axes(vdiff, 1)
        v = vdiff[i, :]
        n = normal[i, :]
        vnorm = norm(v)
        force_dir = maximum(abs.(v)) == 0 ? zeros(FT, 2) : v/vnorm
        friction = G * Δl[i] * Δt * vnorm * -dot(force_dir, v) * force_dir
        if norm(friction) > consts.μ*norm(n)
            friction = -consts.μ*norm(n)*force_dir
        end
        force[i, :] = friction
    end
    return force
end

"""
    floe_floe_interaction!(
        ifloe,
        i,
        jfloe,
        j,
        nfloes,
        consts,
        Δt,
    )

If the two floes interact, update floe i's interactions accordingly. Floe j is
not update here so that the function can be parallelized.
Inputs:
    ifloe       <Floe> first floe in potential interaction
    i           <Int> index of ifloe in model's list of floes 
    jfloe       <Floe> second floe in potential interaction
    j           <Int> index of jfloe in model's list of floes 
    nfloes      <Int> number of non-ghost floes in the simulation this timestep
    consts      <Constants> model constants needed for calculations
    Δt          <Int> Simulation's current timestep
    max_overlap <Float> Percent two floes can overlap before marking them
                        for combination with remove/transfer
Outputs:
    remove   <Int> index of floe to remove from simulation due to overlap (will
                be transfered to over floe in collision) - if 0 then no floes to
                be removed from this collision
    transfer <Int> index of floe to transfer mass to if the other floe is to be
                removed - if 0 then no floes transfered from this collision
Note:
    If ifloe interacts with jfloe, only ifloe's interactions field is updated
    with the details of each region of overlap. The interactions field will have
    the following form for each region of overlap with the boundary:
    [Inf, xforce, yforce, xfpoints, yfpoints, overlaps] where the xforce and
    yforce are the forces, xfpoints and yfpoints are the location of the force
    and overlaps is the overlap between the floe and boundary. The overlaps
    field is also added to the floe's overarea field that describes the total
    overlapping area at any timestep. 
"""
function floe_floe_interaction!(
    ifloe,
    i,
    jfloe,
    j,
    nfloes,
    consts,
    Δt,
    max_overlap::FT,
) where {FT<:AbstractFloat}
    remove = Int(0)
    transfer = Int(0)
    ifloe_poly = LG.Polygon(ifloe.coords)
    jfloe_poly = LG.Polygon(jfloe.coords)
    inter_floe = LG.intersection(ifloe_poly, jfloe_poly)
    inter_regions = LG.Polygon[]
    region_areas = FT[]
    total_area = FT(0)
    for geom in LG.getGeometries(inter_floe)
        a = LG.area(geom)::FT
        if a > 0
            push!(inter_regions, geom)
            push!(region_areas, a)
            total_area += a
        end
    end
    if total_area > 0
        # Floes overlap too much - remove floe or transfer floe mass
        if total_area/ifloe.area > max_overlap
            if i <= nfloes  # If i is not a ghost floe
                remove = i
                transfer = j  # Will transfer mass to jfloe
            elseif j <= nfloes  # If j is not a ghost floe
                # Will transfer mass to ifloe - can't be updated here due to parallelization
                remove = j 
            end
        elseif total_area/jfloe.area > max_overlap
            if j <= nfloes  # If j is not a ghost floe
                # Will transfer mass to ifloe - can't be updated here due to parallelization
                remove = j
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
            normal_forces, fpoints, overlaps, Δl = calc_elastic_forces(
                ifloe.coords,
                jfloe.coords,
                inter_regions,
                region_areas,
                force_factor
            )
            #= Calculate frictional forces at each force point - based on
            velocities at force points =#
            np = size(fpoints, 1)
            if np > 0
                iv = repeat([ifloe.u ifloe.v], outer = np) .+
                    ifloe.ξ*(fpoints .- repeat(ifloe.centroid', outer = np)) 
                jv = repeat([jfloe.u jfloe.v], outer = np) .+
                    jfloe.ξ*(fpoints .- repeat(jfloe.centroid', outer = np))
                friction_forces = calc_friction_forces(
                    iv,
                    jv,
                    normal_forces,
                    Δl,
                    consts,
                    Δt
                )
                # Calculate total forces and update ifloe's interactions
                forces = normal_forces .+ friction_forces
                if sum(abs.(forces)) != 0
                    ifloe.interactions = [ifloe.interactions;
                        fill(j, np) forces fpoints zeros(np) overlaps]
                    ifloe.overarea += sum(overlaps)
                end
            end
        end
    end
    return remove, transfer
end

"""
    floe_domain_element_interaction!(floe, boundary, _, _,)

If given floe insersects with an open boundary, the floe is set to be removed
from the simulation.
Inputs:
    floe            <Floe> floe interacting with boundary
    boundary        <OpenBoundary> coordinates of boundary
    _               <Constants> model constants needed in other methods of this
                        function - not needed here
    _               <Int> current simulation timestep - not needed here
    -               <Float> maximum overlap between floe and domain elements -
                        not needed here
Output:
    None. If floe is interacting with the boundary, floe.alive field is set to
    0. Else, nothing is changed. 
"""
function floe_domain_element_interaction!(
    floe,
    boundary::OpenBoundary,
    consts,
    Δt,
    max_overlap,
)
    floe_poly = LG.Polygon(floe.coords)
    bounds_poly = LG.Polygon(boundary.coords)
    # Check if the floe and boundary actually overlap
    if LG.area(LG.intersection(floe_poly, bounds_poly)) > 0
        floe.alive = false
    end
    return
end

"""
    floe_domain_element_interaction!(
        floe,
        ::PeriodicBoundary,
        consts,
        Δt,
    )

If a given floe intersects with a periodic boundary, nothing happens at this
point. Periodic floes pass through boundaries using ghost floes.
Inputs:
        None are used. 
Output:
        None. This function does not do anyting. 
"""
function floe_domain_element_interaction!(
    floe,
    ::PeriodicBoundary,
    consts,
    Δt,
    max_overlap,
)
    return
end

"""
    normal_direction_correct!(
        forces,
        fpoints,
        boundary::AbstractBoundary{North, <:AbstractFloat},
    )

Zero-out forces that point in direction not perpendicular to North boundary wall.
Inputs:
    force       <Array{Float, n, 2}> normal forces on each of the n regions
                    greater than a minimum area
    fpoint      <Array{Float, n, 2}> point force is applied on each of the n
                    regions greater than a minimum area
    boundary    <AbstractBoundary{North, <:AbstractFloat}> domain's northern
                    boundary
Outputs:
    None. All forces in the x direction set to 0 if the point the force is
    applied to is in the northern boundary.
"""
function normal_direction_correct!(
    forces::Matrix{FT},
    fpoints,
    boundary::AbstractBoundary{North, <:AbstractFloat},
) where {FT}
    forces[fpoints[:, 2] .>= boundary.val, 1] .= FT(0.0)
    return
end

"""
    normal_direction_correct!(
        forces,
        fpoints,
        boundary::AbstractBoundary{South, <:AbstractFloat},
    )

Zero-out forces that point in direction not perpendicular to South boundary wall.
See normal_direction_correct! on northern wall for more information
"""
function normal_direction_correct!(
    forces::Matrix{FT},
    fpoints,
    boundary::AbstractBoundary{South, <:AbstractFloat},
) where {FT}
        forces[fpoints[:, 2] .<= boundary.val, 1] .= FT(0.0)
        return
    end

"""
    normal_direction_correct!(
        forces,
        fpoints,
        boundary::AbstractBoundary{East, <:AbstractFloat},
    )

Zero-out forces that point in direction not perpendicular to East boundary wall.
See normal_direction_correct! on northern wall for more information
"""
function normal_direction_correct!(
    forces::Matrix{FT},
    fpoints,
    boundary::AbstractBoundary{East, <:AbstractFloat},
) where {FT}
    forces[fpoints[:, 1] .>= boundary.val, 2] .= FT(0.0)
    return
end

"""
    normal_direction_correct!(
        forces,
        fpoints,
        boundary::AbstractBoundary{West, <:AbstractFloat},
    )

Zero-out forces that point in direction not perpendicular to West boundary wall.
See normal_direction_correct! on northern wall for more information
"""
function normal_direction_correct!(
    forces::Matrix{FT},
    fpoints,
    boundary::AbstractBoundary{West, <:AbstractFloat},
) where {FT}
    forces[fpoints[:, 1] .<= boundary.val, 2] .= FT(0.0)
    return
end

"""
    normal_direction_correct!(
        forces,
        fpoints,
        ::TopographyElement,
    )

No forces should be zero-ed out in collidions with topography elements. 
Inputs:
        None used.
Outputs:
        None.
"""
function normal_direction_correct!(
    forces,
    fpoints,
    ::TopographyElement,
)
    return
end

"""
    floe_domain_element_interaction!(
        floe,
        element,
        consts,
        Δt,
    )

If floe intersects with given element (either collision boundary or
topography element), floe interactions field and overarea field are updated.
Inputs:
    floe            <Floe> floe interacting with element
    element         <Union{CollisionBoundary, TopographyElement}> coordinates of
                        element
    consts          <Constants> model constants needed for calculations
    Δt              <Int> current simulation timestep
    max_overlap     <Float> Percent a floe can overlap with a collision wall
                            or topography before being killed/removed
Outputs:
    None. If floe interacts, the floe's interactions field is updated with the
    details of each region of overlap. The interactions field will have the
    following form for each region of overlap with the element:
    [Inf, xforce, yforce, xfpoints, yfpoints, overlaps] where the xforce and
    yforce are the forces, xfpoints and yfpoints are the location of the force
    and overlaps is the overlap between the floe and element. The overlaps field
    is also added to the floe's overarea field that describes the total
    overlapping area at any timestep.
"""
function floe_domain_element_interaction!(
    floe,
    element::Union{CollisionBoundary, TopographyElement},
    consts,
    Δt,
    max_overlap::FT,
) where {FT}
    floe_poly = LG.Polygon(floe.coords)
    bounds_poly = LG.Polygon(element.coords)
    # Check if the floe and element actually overlap
    inter_floe = LG.intersection(floe_poly, bounds_poly)
    inter_regions = LG.Polygon[]
    region_areas = FT[]
    max_area = FT(0)
    for geom in LG.getGeometries(inter_floe)
        a = LG.area(geom)::FT
        if a > 0
            push!(inter_regions, geom)
            push!(region_areas, a)
            if a > max_area
                max_area = a
            end
        end
    end
    if max_area > 0
        # Regions overlap too much
        if max_area/floe.area > max_overlap
            floe.alive = false
        else
            # Constant needed for force calculations
            force_factor = consts.E * floe.height / sqrt(floe.area)
            # Calculate normal forces, force points, and overlap areas
            normal_forces, fpoints, overlaps, Δl =  calc_elastic_forces(
                floe.coords,
                element.coords,
                inter_regions,
                region_areas,
                force_factor,
            )
            normal_direction_correct!(normal_forces, fpoints, element)
            # Calculate frictional forces at each force point
            np = size(fpoints, 1)
            vfloe = repeat([floe.u floe.v], outer = np) .+
                floe.ξ*(fpoints .- repeat(floe.centroid', outer = np)) 
            vbound = repeat(zeros(FT, 1, 2), outer = np)
            friction_forces = calc_friction_forces(
                vfloe,
                vbound,
                normal_forces,
                Δl,
                consts,
                Δt,
            )
            # Calculate total forces and update ifloe's interactions
            forces = normal_forces .+ friction_forces
            if sum(abs.(forces)) != 0
                floe.interactions = [floe.interactions;
                    fill(Inf, np) forces fpoints zeros(np) overlaps]
                floe.overarea += sum(overlaps)
            end
        end
    end
    return
end

"""
    floe_domain_interaction!(
        floe,
        domain::DT,
        consts,
        max_overlap,
    )

If the floe interacts with the domain, update the floe accordingly. Dispatches
on different boundary types within the domain.
Inputs:
    floe        <Floe> floe interacting with boundary
    domain      <Domain> model domain
    consts      <Constants> model constants needed for calculations
    Δt          <Int> current simulation timestep
    max_overlap <Float> Percent a floe can overlap with a collision wall
                        or topography before being killed/removed
Outputs:
    None. Floe is updated according to which boundaries it interacts with and
    the types of those boundaries. 
"""
function floe_domain_interaction!(
    floe,
    domain::Domain,
    consts,
    Δt,
    max_overlap::FT,
) where {FT<:AbstractFloat}
    centroid = floe.centroid
    rmax = floe.rmax
    nbound = domain.north
    sbound = domain.south
    ebound = domain.east
    wbound = domain.west

    if centroid[2] + rmax > nbound.val
        floe_domain_element_interaction!(
            floe,
            nbound,
            consts,
            Δt,
            max_overlap,
        )
    end
    if centroid[2] - rmax < sbound.val
        floe_domain_element_interaction!(
            floe,
            sbound,
            consts,
            Δt,
            max_overlap,
        )
    end
    if centroid[1] + rmax > ebound.val
        floe_domain_element_interaction!(
            floe,
            ebound,
            consts,
            Δt,
            max_overlap,
        )
    end
    if centroid[1] - rmax < wbound.val
        floe_domain_element_interaction!(
            floe,
            wbound,
            consts,
            Δt,
            max_overlap,
        )
    end

    for topo_element in domain.topography
        if sum((topo_element.centroid .- floe.centroid).^2) < (topo_element.rmax + floe.rmax)^2
            floe_domain_element_interaction!(
                floe,
                topo_element,
                consts,
                Δt,
                max_overlap,
            )
        end
    end

    return
end

"""
    calc_torque!(floe)

Calculate a floe's torque based on the interactions.
Inputs:
        floe  <Floe> floe in model
Outputs:
        None. Floe's interactions field is updated with calculated torque.
"""
function calc_torque!(floe::Union{LazyRow{<:Floe{FT}}, Floe{FT}}) where {FT<:AbstractFloat}
    inters = floe.interactions
    if !isempty(inters)
        dir = [inters[:, xpoint] .- floe.centroid[1] inters[:, ypoint] .- floe.centroid[2] zeros(FT, size(inters, 1))]
        frc = [inters[:, xforce] inters[:, yforce] zeros(FT, size(inters, 1))]
        for i in axes(dir, 1)
            idir = vec(dir[i, :])
            ifrc = vec(frc[i, :])
            itorque = cross(idir, ifrc)
            floe.interactions[i, torque] = itorque[3]
        end
    end
end

"""
    timestep_collisions!(
        floes,
        n_init_floes,
        domain,
        remove,
        transfer,
        consts,
        Δt,
        collision_settings,
        spinlock,
    )

Resolves collisions between pairs of floes and calculates the forces and torques
caused by those collisions.
Inputs:
    floes               <StructArray{Floe}> model's list of floes
    n_init_floes        <Int> number of floes without ghost floes
    domain              <Domain> model's domain
    remove              <Vector{Int}> indices of floes to remove from list of
                            floes due to floe-floe overlap
    transfer            <Vector{Int}> indices of floes to transfer mass to from
                            list of floes due to floe-floe overlap
    consts              <Constants> simulation constants
    Δt                  <Int> length of simulation timestep in seconds
    collision_settings  <CollisionSettings> simulation collision settings
    spinlock            <Thread.SpinLock>
"""
function timestep_collisions!(
    floes::StructArray{<:Floe{FT}},
    n_init_floes,
    domain,
    remove,
    transfer,
    consts,
    Δt,
    collision_settings,
    spinlock,
) where {FT<:AbstractFloat}
    collide_pairs = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    # floe-floe collisions for floes i and j where i<j
    for i in eachindex(floes)
        # reset collision values
        fill!(floes.collision_force[i], FT(0))
        floes.collision_trq[i] = FT(0.0)
        floes.interactions[i] = zeros(FT, 0, 7)
        for j in i+1:length(floes)
            id_pair, ghost_id_pair =
                if floes.id[i] > floes.id[j]
                    (floes.id[i], floes.id[j]), (floes.ghost_id[i], floes.ghost_id[j])
                else
                    (floes.id[j], floes.id[i]), (floes.ghost_id[j], floes.ghost_id[i])
                end
            # If checking two distinct floes (i.e. not a parent ghost pair) and they are in close proximity continue
            if (id_pair[1] != id_pair[2]) && (sum((floes.centroid[i] .- floes.centroid[j]).^2) < (floes.rmax[i] + floes.rmax[j])^2)
                # Never seen any combo of these floes/ghosts
                new_collision = !(id_pair in keys(collide_pairs))
                if new_collision
                    Threads.lock(spinlock) do
                        collide_pairs[id_pair] = ghost_id_pair
                    end
                end
                # New collision or floe and ghost colliding with same floe - not a repeat collision
                if new_collision || (ghost_id_pair[1] == collide_pairs[id_pair][1]) ⊻ (ghost_id_pair[2] == collide_pairs[id_pair][2])
                    iremove, itransfer = floe_floe_interaction!(
                        LazyRow(floes, i),
                        i,
                        LazyRow(floes, j),
                        j,
                        n_init_floes,
                        consts,
                        Δt,
                        collision_settings.floe_floe_max_overlap,
                    )
                    if iremove != 0 || itransfer != 0
                        remove[i] = iremove
                        transfer[i] = itransfer
                    end
                end
            end
        end
        floe_domain_interaction!(
            LazyRow(floes, i),
            domain,
            consts,
            Δt,
            collision_settings.floe_domain_max_overlap,
        )
    end
    # Update floes not directly calculated above where i>j - can't be parallelized
    for i in eachindex(floes)
        # if we are removing floe j = remove[i] floe j's mass must be transfered to floe i
        if i <= n_init_floes && remove[i] > 0 && remove[i] != i
            transfer[remove[i]] = i
        end
        ij_inters = floes.interactions[i]
        if !isempty(ij_inters)
            for inter_idx in axes(ij_inters, 1)  # Loop over each interaction with Floe i
                j = ij_inters[inter_idx, floeidx]  # Index of floe to update in model floe list
                if j <= length(floes) && j > i
                    jidx = Int(j)
                    floes.interactions[jidx] = [floes.interactions[jidx]; i -ij_inters[inter_idx, xforce] -ij_inters[inter_idx, yforce] #=
                                             =# ij_inters[inter_idx, xpoint] ij_inters[inter_idx, ypoint] FT(0.0) ij_inters[inter_idx, overlap]]
                    floes.overarea[jidx] += ij_inters[inter_idx, overlap]
                end
            end
        end
    end
    # Calculate total force and torque on parents by summing interactins on parent and children
    for i in range(1, n_init_floes)
        for g in floes.ghosts[i]
            g_inters = floes.interactions[g]
            g_inters[:, xpoint] .-= (floes.centroid[g][1] - floes.centroid[i][1])
            g_inters[:, ypoint] .-= (floes.centroid[g][2] - floes.centroid[i][2])
            floes.interactions[i] = [floes.interactions[i]; g_inters]
        end
        calc_torque!(LazyRow(floes, i))
        floes.collision_force[i][1] += sum(@view floes.interactions[i][:, xforce])
        floes.collision_force[i][2] += sum(@view floes.interactions[i][:, yforce])
        floes.collision_trq[i] += sum(@view floes.interactions[i][:, torque])
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
    if LG.area(LG.intersection(LG.Polygon(element.coords), LG.Polygon(boundary.coords))) > 0
        # ghosts of existing ghosts and original element
        append!(new_ghosts, deepcopy.(ghosts))
        push!(new_ghosts, deepcopy(element))
        for i in eachindex(new_ghosts)
            translate!(new_ghosts.coords[i], trans_vec[1], trans_vec[2])
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
        if (elem.centroid[1] - elem.rmax < wbound.val)  # passing through western boundary
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
    add_floe_ghosts!(floes, max_boundary, min_boundary)

Add ghosts of all of the given floes passing through the two given boundaries to the list of floes.
Inputs:
        floes           <StructArray{Floe{FT}}> list of floes to find ghosts for
        max_boundary    <PeriodicBoundary> either northern or eastern boundary  of domain
        min_boundary    <PeriodicBoundary> either southern or western boundary of domain
Outputs:
        None. Ghosts of floes are added to floe list. 
"""
function add_floe_ghosts!(floes::StructArray{Floe{FT}}, max_boundary, min_boundary) where {FT <: AbstractFloat}
    nfloes = length(floes)
    for i in eachindex(floes)  # uses initial length of floes so we can append to list
        f = floes[i]
        if f.alive && f.ghost_id == 0  # the floe is alive and a parent floe
            f, new_ghosts = find_ghosts(f, floes[f.ghosts], max_boundary, min_boundary, FT)
            if !isempty(new_ghosts)
                nghosts = length(new_ghosts)
                new_ghosts.ghost_id .= range(1, nghosts) .+ length(f.ghosts)
                empty!.(new_ghosts.ghosts)  # remove ghost floes ghosts as these were added to parent
                append!(floes, new_ghosts)  # add ghosts to floe list
                append!(f.ghosts, (nfloes + 1):(nfloes + nghosts))  # index of ghosts floes saved in parent
                nfloes += nghosts
                floes[i] = f
            end
        end
    end
    return
end

"""
    add_ghosts!(
        elems,
        domain,
    )

When there are no periodic boundaries, no ghosts should be added.
Inputs:
        None are used. 
Outputs:
        None. 
"""
function add_ghosts!(
    elems,
    ::Domain{
        FT,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
    },
) where {FT<:AbstractFloat}
    return
end

"""
    add_ghosts!(
        elems,
        domain,
    )

Add ghosts for elements that pass through the northern or southern boundaries.
Inputs:
        elems   <StructArray{Floe} or StructArray{TopographyElement}> list of
                    elements to add ghosts to
        domain  <Domain{
                    Float,
                    PeriodicBoundary,
                    PeriodicBoundary,
                    NonPeriodicBoundary,
                    NonPeriodicBoundary,
                }> domain with northern and southern periodic boundaries
Outputs:
        None. Ghosts are added to list of elements.
"""
function add_ghosts!(
    elems,
    domain::Domain{
        FT,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
    },
) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.north, domain.south)
    return
end

"""
    add_ghosts!(
        elems,
        domain,
    )

Add ghosts for elements that pass through the eastern or western boundaries. 
Inputs:
    elems   <StructArray{Floe} or StructArray{TopographyElement}> list of
                elements to add ghosts to
    domain  <Domain{
                Float,
                NonPeriodicBoundary,
                NonPeriodicBoundary,
                PeriodicBoundary,
                PeriodicBoundary,
            }> domain with eastern and western periodic boundaries 
Outputs:
    None. Ghosts are added to list of elements.
"""
function add_ghosts!(
    elems,
    domain::Domain{
        FT,
        <:NonPeriodicBoundary,
        <:NonPeriodicBoundary,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
    },
) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.east, domain.west)
    return
end

"""
    add_ghosts!(
        elems,
        domain,
    )

Add ghosts for elements that pass through any of the boundaries. 
Inputs:
    elems   <StructArray{Floe} or StructArray{TopographyElement}> list of
                elements to add ghosts to
    domain  <Domain{
                AbstractFloat,
                PeriodicBoundary,
                PeriodicBoundary,
                PeriodicBoundary,
                PeriodicBoundary,
            }> domain with all boundaries
Outputs:
        None. Ghosts are added to list of elements.
"""
function add_ghosts!(
    elems,
    domain::Domain{
        FT,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
        <:PeriodicBoundary,
        <:PeriodicBoundary
    }
) where {FT<:AbstractFloat}
    add_floe_ghosts!(elems, domain.east, domain.west)
    add_floe_ghosts!(elems, domain.north, domain.south)
    return
end
