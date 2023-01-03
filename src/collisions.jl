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
        mag = sqrt(Δx^2 + Δy^2)
        Δl = mag
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

        # Need to find which new overlap region corresponds to region k
        new_region_overlaps = LG.getGeometries(LG.intersection(new_inter_floe, region))
        if !isempty(new_region_overlaps)
            ~, max_idx = findmax(LG.area, new_region_overlaps)
            new_overlap_area = LG.area(LG.getGeometries(new_inter_floe)[max_idx])
            if new_overlap_area/area > 1  # Force increased overlap
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
        consts          <Constants> model constants needed for calculations
        t               <Type> Float type model is running on (Float64 or Float32)
Outputs:
        force   <Array{Float, n, 2}> normal forces on each of the n regions greater than a minimum area
        fpoint  <Array{Float, n, 2}> point force is applied on each of the n regions greater than a minimum area
        overlap <Array{Float, n, 2}> area of each of the n regions greater than a minimum area
        Δl      <Float> mean length of distance between intersection points
"""
function calc_elastic_forces(c1, c2, regions, region_areas, force_factor, consts, t::Type{T} = Float64) where T
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
        Δl = T(0.0)
        for k in 1:ncontact
            normal_force = zeros(T, 1, 2)
            if region_areas[k] != 0
                cx, cy = LG.GeoInterface.coordinates(LG.centroid(regions[k]))::Vector{Float64}
                fpoint[k, :] = [cx, cy]
                normal_force, Δl = calc_normal_force(c1, c2, regions[k], region_areas[k], ipoints, force_factor, T)
            end
            force[k, :] = normal_force
        end
        return force, fpoint, overlap, Δl
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
        vnorm = norm(v)
        force_dir = maximum(abs.(v)) == 0 ? zeros(T, 2) : v/vnorm
        friction = G * Δl * Δt * vnorm * -dot(force_dir, v) * force_dir
        if norm(friction) > consts.μ*norm(normal)
            friction = consts.μ*norm(normal)*force_dir
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
            if i <= nfloes
                remove = i
                transfer = j
            elseif j <= nfloes
                remove = j  # Will transfer mass to ifloe - can't be updated here due to parallelization
            end
        elseif total_area/jfloe.area > 0.55
            if j <= nfloes
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
                                        inter_regions, region_areas, force_factor, consts, T)
            # Calculate frictional forces at each force point - based on velocities at force points
            np = size(fpoints, 1)
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
    return remove, transfer
end

"""
    floe_domain_element_interaction!(floe, boundary, _, _, _)

If given floe insersects with an open boundary, the floe is set to be removed from the simulation.
Inputs:
        floe            <Floe> floe interacting with boundary
        boundary        <OpenBoundary> coordinates of boundary
        _               <Constants> model constants needed in other methods of this function - not needed here
        _               <Int> current simulation timestep 
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
    normal_direction_correct!(forces, fpoints, boundary::AbstractBoundary{North, <:AbstractFloat}, ::Type{T} = Float64)

Zero-out forces that point in direction not perpendicular to North or South boundary wall.
Inputs:
        force       <Array{Float, n, 2}> normal forces on each of the n regions greater than a minimum area
        fpoint      <Array{Float, n, 2}> point force is applied on each of the n regions greater than a minimum area
        boundary    <AbstractBoundary{North, <:AbstractFloat}> domain's northern boundary 
                    <Type> Float type model is running on (Float64 or Float32) - not needed here
Outputs: None. All forces in the x direction set to 0 if the point the force is applied is the northern or southern boundary value.
"""
function normal_direction_correct!(forces, fpoints, boundary::Union{AbstractBoundary{North, <:AbstractFloat},
AbstractBoundary{South, <:AbstractFloat}}, ::Type{T} = Float64) where T
    forces[fpoints[:, 2] .== boundary.val, 1] .= T(0.0)
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
function normal_direction_correct!(forces, fpoints, boundary::Union{AbstractBoundary{East, <:AbstractFloat},
                                                                    AbstractBoundary{West, <:AbstractFloat}}, ::Type{T} = Float64) where T
    forces[fpoints[:, 1] .== boundary.val, 2] .= T(0.0)
    return
end

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
            return zeros(T, 1, 2), zeros(T, 1, 2), fill(T(Inf), 1)
        end
        # Constant needed for force calculations
        force_factor = consts.E * floe.height / sqrt(floe.area)
        # Calculate normal forces, force points, and overlap areas
        normal_forces, fpoints, overlaps, Δl =  calc_elastic_forces(floe.coords, element.coords,
                                        inter_regions, region_areas, force_factor, consts, T)
        # Calculate frictional forces at each force point - based on velocities at force points
        np = size(fpoints, 1)
        vfloe = repeat([floe.u floe.v], outer = np) .+ floe.ξ*(fpoints .- repeat(floe.centroid', outer = np)) 
        vbound = repeat(zeros(T, 1, 2), outer = np)
        friction_forces = calc_friction_forces(vfloe, vbound, normal_forces, Δl, consts, Δt, T)
        # Calculate total forces and update ifloe's interactions
        forces = normal_forces .+ friction_forces
        if sum(abs.(forces)) != 0
            normal_direction_correct!(forces, fpoints, element, T)
            floe.interactions = [floe.interactions; fill(Inf, np) forces fpoints zeros(np) overlaps']
            floe.overarea += sum(overlaps)
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

    if centroid[1] > ebound.val || centroid[1] < wbound.val || centroid[2] > nbound.val || centroid[2] < sbound.val
        floe.alive = 0
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
    if size(inters, 1) > 1
        dir = [inters[:, "xpoint"] .- floe.centroid[1] inters[:, "ypoint"] .- floe.centroid[2] zeros(T, size(inters, 1))]
        frc = [inters[:, "xforce"] inters[:, "yforce"] zeros(T, size(inters, 1))]
        for i in axes(dir, 1)
            idir = vec(dir[i, :])
            ifrc = vec(frc[i, :])
            itorque = cross(idir, ifrc)
            floe.interactions[i, "torque"] = itorque[3]
        end
    end
end

