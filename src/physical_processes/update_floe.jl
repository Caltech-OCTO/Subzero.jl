"""
Functions to update floe's shape and fields. Function typically called from
other physical processes, other than timestep_floe_properties!, which is called
from the timestep_simulation! function.
"""

"""
    replace_floe!(
        floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
        new_poly,
        new_mass,
        consts,
        npoints,
        rng,
    )
Updates existing floe shape and related physical properties based of the polygon
defining the floe.
Inputs:
    floe        <Union{Floe, LazyRow{Floe}}> floe to update
    new_poly    <LG.Polygon> polygon representing new outline of floe
    new_mass    <AbstractFloat> mass of floe
    consts      <Constants> simulation's constants
    npoints        <Int> number of monte carlo points to attempt to generate
    rng         <RNG> random number generator
Ouputs:
    Updates a given floe's physical properties given new shape and total mass.
"""
function replace_floe!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    new_poly,
    new_mass,
    consts,
    coupling_settings,
    Δg,
    rng,
) where {FT}
    # Floe shape
    floe.centroid = find_poly_centroid(new_poly)
    floe.coords = find_poly_coords(new_poly)::PolyVec{FT}
    floe.coords = [orient_coords(floe.coords[1])]
    if floe.coords[1][1] != floe.coords[1][end]
        push!(floe.coords, floe.coords[1][1])
    end
    floe.area = LG.area(new_poly)
    floe.height = new_mass/(floe.area * consts.ρi)
    floe.mass = new_mass
    floe.moment = calc_moment_inertia(
        floe.coords,
        floe.centroid,
        floe.height;
        ρi = consts.ρi,
    )
    floe.angles = calc_poly_angles(floe.coords)
    floe.α = FT(0)
    translate!(floe.coords, -floe.centroid[1], -floe.centroid[2])
    floe.rmax = sqrt(maximum([sum(c.^2) for c in floe.coords[1]]))
    # Floe monte carlo points
    x_subpoints, y_subpoints, status = generate_floe_points(
        FT,
        coords,
        rmax,
        area,
        status,
        Δg,
        coupling_settings,
        rng,
    )

    translate!(floe.coords, floe.centroid[1], floe.centroid[2])
    floe.x_subpoints = x_subpoints
    floe.y_subpoints = y_subpoints
    # Floe status / identification
    floe.status = status
    return
end

"""
    conserve_momentum_combination!(
        mass_tmp,
        moment_tmp,
        x_tmp,
        y_tmp,
        Δt,
        keep_floe,
        remove_floe = nothing,
    )
Update current and previous velocity/acceleration fields to conserve momentum of
a floe whose shape has been changed, given the previous mass, momentum, and
centroid.
Inputs:
    mass_tmp    <AbstractFloat> original mass of floe before shape change
    moment_tmp  <AbstractFloat> original moment of intertia of floe before shape
                    change
    x_tmp       <AbstractFloat> original x-coordinate of centroid of floe before
                    shape change
    y_tmp       <AbstractFloat> original y-coordinate of centroid of floe before
                    shape change
    Δt          <Int> timestep of simulation in seconds
    keep_floe   <Union{Floe, LazyRow{Floe}}> floe whose shape has been changed
    remove_floe <Union{Floe, LazyRow{Floe}}> if keep_floe's shape has been
                    changed due to an interaction with another floe, remove_floe
                    is that floe - optional parameter
Output:
    None. keep_floe's u, v, ξ, p_dxdt, p_dydt, p_dαdt, p_dudt, p_dvdt, and
    p_dξdt fields all updated to preserve momentum. 
Note:
    Function depends on conservation of mass. That is, if we do not have
    remove_floe then mass_tmp must be equal to keep_floe, as the mass has not
    changed, or if we have remove_floe then the sum of remove_floe's mass and
    mass_tmp must be equal to keep_floe's mass.
"""
function conserve_momentum_combination!(
    mass_tmp,
    moment_tmp,
    x_tmp,
    y_tmp,
    Δt,
    keep_floe,
    remove_floe = nothing,
)
    # Contributions from first floe in current timestep
    keep_floe.ξ = keep_floe.ξ * moment_tmp +  # change in spin
        mass_tmp * (  # change in orbital
            keep_floe.v * (x_tmp - keep_floe.centroid[1]) +
            keep_floe.u * (keep_floe.centroid[2] - y_tmp)
        )
    keep_floe.u = keep_floe.u * mass_tmp
    keep_floe.v = keep_floe.v * mass_tmp
    # Contributions from first floe in previous timestep
    keep_floe.p_dαdt = keep_floe.p_dαdt * moment_tmp +  # change in spin
        mass_tmp * (  # change in orbital
            keep_floe.p_dydt * (x_tmp - keep_floe.centroid[1]) +
            keep_floe.p_dxdt * (keep_floe.centroid[2] - y_tmp)
        )
    keep_floe.p_dxdt = keep_floe.p_dxdt * mass_tmp
    keep_floe.p_dydt = keep_floe.p_dydt * mass_tmp

    # Contributions from any other floes combining with first
    if !isnothing(remove_floe)
        # Current timestep
        keep_floe.ξ += remove_floe.ξ * remove_floe.moment +
            remove_floe.mass * (
                remove_floe.v * (remove_floe.centroid[1] - keep_floe.centroid[1]) +
                remove_floe.u * (keep_floe.centroid[2] - remove_floe.centroid[2])
            )
        keep_floe.u += remove_floe.u * remove_floe.mass
        keep_floe.v += remove_floe.v * remove_floe.mass
        # Previous timestep
        keep_floe.p_dαdt += remove_floe.p_dαdt * remove_floe.moment +
            remove_floe.mass * (
                remove_floe.p_dydt * (remove_floe.centroid[1] - keep_floe.centroid[1]) +
                remove_floe.p_dxdt * (keep_floe.centroid[2] - remove_floe.centroid[2])
            )
        keep_floe.p_dxdt += remove_floe.p_dxdt * remove_floe.mass
        keep_floe.p_dydt += remove_floe.p_dydt * remove_floe.mass
    end

    # Average added values
    keep_floe.ξ /= keep_floe.moment
    keep_floe.u /= keep_floe.mass
    keep_floe.v /= keep_floe.mass
    keep_floe.p_dαdt /= keep_floe.moment
    keep_floe.p_dxdt /= keep_floe.mass
    keep_floe.p_dydt /= keep_floe.mass
    # Calculate previous accelerations
    keep_floe.p_dξdt = (keep_floe.ξ - keep_floe.p_dαdt) / Δt
    keep_floe.p_dudt = (keep_floe.u - keep_floe.p_dxdt) / Δt
    keep_floe.p_dvdt = (keep_floe.v - keep_floe.p_dydt) / Δt
    return
end

"""
    conserve_momentum_fracture!(
        init_floe,
        new_floes,
        Δt,
    )
Update new_floes's current and previous velocity/acceleration fields to conserve
momentum when a floe has been fractured into several new floes, given the
previous mass, momentum, and centroid. The assumption is made that each new floe
has the same velocities/accelerations
Inputs:
    init_floe   <Union{Floe, LazyRow{Floe}}> original floe
    new_floes   <StructArray{Floe}> fractured pieces of original floe
    Δt          <Int> simulation's timestep in seconds
Output:
    None. new_floes velocities and accelerations are updated for current and
    previous timestep to conserve momentum.
"""
function conserve_momentum_fracture!(
    init_floe,
    new_floes,
    Δt,
)
    if !isempty(new_floes)
        x_init, y_init = init_floe.centroid
        # conserve linear momentum by keeping linear velocities the same
        new_floes.u .= init_floe.u 
        new_floes.v .= init_floe.v
        new_floes.p_dxdt .= init_floe.p_dxdt
        new_floes.p_dydt .= init_floe.p_dydt
        new_floes.p_dudt .= init_floe.p_dudt
        new_floes.p_dvdt .= init_floe.p_dvdt
        # conserve rotational moment by offsetting change in orbital momentum
        sum_moments = sum(new_floes.moment)
        new_floes.ξ .= init_floe.moment * init_floe.ξ + # initial spin velocity
            init_floe.mass * ( # initial orbital velocity
                x_init * init_floe.v -
                y_init * init_floe.u
            )
        new_floes.p_dαdt .= init_floe.moment * init_floe.p_dαdt + 
            init_floe.mass * (
                x_init * init_floe.p_dydt -
                y_init * init_floe.p_dxdt
            )
        for i in eachindex(new_floes)
            new_floes.ξ .-= new_floes.mass[i] * (
                new_floes.centroid[i][1] * new_floes.v[i] -
                new_floes.centroid[i][2] * new_floes.u[i]
            )
            new_floes.p_dαdt .-= new_floes.mass[i] * (
                new_floes.centroid[i][1] * new_floes.p_dydt[i] -
                new_floes.centroid[i][2] * new_floes.p_dxdt[i]
            )
        end
        new_floes.ξ ./= sum_moments
        new_floes.p_dαdt ./= sum_moments
        # Calculate previous rotational acceleration
        new_floes.p_dξdt .= (new_floes.ξ[1] - new_floes.p_dαdt[1]) / Δt
    end
end

"""
    calc_stress!(floe)

Calculates the stress on a floe for current collisions given interactions and
floe properties.
Inputs:
    inters      <Matrix{AbstractFloat}> matrix of floe interactions
    centroid    <Vector{AbstractFloat}> floe centroid as [x, y] coordinates
    area        <AbstractFloat> floe area
    height      <AbstractFloat> floe height
Outputs:
    Caculates stress on floe at current timestep from interactions
"""
function calc_stress!(floe)
    # Stress calcultions
    xi, yi = floe.centroid
    inters = floe.interactions
    # Calculates timestep stress
    stress = fill(1/(2* floe.area * floe.height), 2, 2)
    stress[1, 1] *= sum((inters[:, xpoint] .- xi) .* inters[:, xforce]) +
        sum(inters[:, xforce] .* (inters[:, xpoint] .- xi))
    stress[1, 2] *= sum((inters[:, ypoint] .- yi) .* inters[:, xforce]) + 
        sum(inters[:, yforce] .* (inters[:, xpoint] .- xi))
    stress[2, 1] *= sum((inters[:, xpoint] .- xi) .* inters[:, yforce]) +
        sum(inters[:, xforce] .* (inters[:, ypoint] .- yi))
    stress[2, 2] *= sum((inters[:, ypoint] .- yi) .* inters[:, yforce]) +
        sum(inters[:, yforce] .* (inters[:, ypoint] .- yi))
    # Add timestep stress to stress history
    push!(floe.stress_history, stress)
    # Average stress history to find floe's average stress
    floe.stress = mean(floe.stress_history)
    return
end

"""
    calc_strain!(coords, centroid, u, v, ξ, area)

Calculates the strain on a floe given the velocity at each vertex
Inputs:
    floe        <Floe{AbstractFloat}> a floe
Outputs:
    strain      <Matrix{AbstractFloat}> 2x2 matrix for floe strain 
"""
function calc_strain!(floe)
    # coordinates of floe centered at centroid
    translate!(floe.coords, -floe.centroid[1], -floe.centroid[2])
    xcoords, ycoords = separate_xy(floe.coords)
    translate!(floe.coords, floe.centroid[1], floe.centroid[2])
    # Find distance between each vertex
    if xcoords[1] != xcoords[end] && ycoords[1] != ycoords[end]
        push!(xcoords, xcoords[1])
        push!(ycoords, ycoords[1])
    end
    xcoords_diff = diff(xcoords)
    ycoords_diff = diff(ycoords)
    # u and v velocities of floes at each vertex
    ucoords = fill(floe.u, size(xcoords))
    vcoords = fill(floe.v, size(xcoords))
    for i in eachindex(ucoords)
        rad = sqrt(xcoords[i]^2 + ycoords[i]^2)
        θ = atan(ycoords[i], xcoords[i])
        ucoords[i] -= floe.ξ * rad * sin(θ)
        vcoords[i] += floe.ξ * rad * cos(θ)
    end
    ucoords_diff = diff(ucoords)
    vcoords_diff = diff(vcoords)
    fill!(floe.strain, 1/(2 * floe.area))
    floe.strain[1, 1] *= sum(ucoords_diff .* ycoords_diff) # dudx
    floe.strain[1, 2] *= 0.5(sum(ucoords_diff .* xcoords_diff) +
        sum(vcoords_diff .* ycoords_diff)) # dudy + dvdx
    floe.strain[2, 1] = floe.strain[1, 2]
    floe.strain[2, 2] *= sum(vcoords_diff .* xcoords_diff) # dvdy
    return
end

"""
    timestep_floe(floe)

Update floe position and velocities using second-order time stepping with
tendencies calculated at previous timesteps. Height, mass, stress, and strain
also updated based on previous timestep thermodynamics and interactions with
other floes. 
Input:
        floe        <Floe>
        Δt          <Int> simulation timestep in second
        max_height  <AbstractFloat> maximum floe height
Output:
        None. Floe's fields are updated with values.
"""
function timestep_floe_properties!(
    floes,
    tstep,
    Δt,
    max_height,
)
    Threads.@threads for i in eachindex(floes)
        cforce = floes.collision_force[i]
        ctrq = floes.collision_trq[i]
        # Update stress
        if !isempty(floes.interactions[i])
            calc_stress!(LazyRow(floes, i))
        end
        # Ensure no extreem values due to model instability
        if floes.height[i] > max_height
            @warn "Reducing height to 10 m"
            floes.height[i] = max_height
        end

        while maximum(abs.(cforce)) > floes.mass[i]/(5Δt)
            @warn "Decreasing collision forces by a factor of 10"
            cforce = cforce ./ 10
            ctrq = ctrq ./ 10
        end
        
        # Update floe based on thermodynamic growth
        h = floes.height[i]
        Δh = floes.hflx_factor[i] / h
        hfrac = (h + Δh) / h
        floes.mass[i] *= hfrac
        floes.moment[i] *= hfrac
        floes.height[i] -= Δh
        h = floes.height[i]

        # Update ice coordinates with velocities and rotation
        Δx = 1.5Δt*floes.u[i] - 0.5Δt*floes.p_dxdt[i]
        Δy = 1.5Δt*floes.v[i] - 0.5Δt*floes.p_dydt[i]
        Δα = 1.5Δt*floes.ξ[i] - 0.5Δt*floes.p_dαdt[i]
        floes.α[i] += Δα

        translate!(
            floes.coords[i],
            -floes.centroid[i][1],
            -floes.centroid[i][2],
        )
        rotate_radians!(floes.coords[i], Δα)
        floes.centroid[i] .+= [Δx, Δy]
        translate!(
            floes.coords[i],
            floes.centroid[i][1],
            floes.centroid[i][2],
        )
        floes.p_dxdt[i] = floes.u[i]
        floes.p_dydt[i] = floes.v[i]
        floes.p_dαdt[i] = floes.ξ[i]

        # Update ice velocities with forces and torques
        dudt = (floes.fxOA[i] + cforce[1])/floes.mass[i]
        dvdt = (floes.fyOA[i] + cforce[2])/floes.mass[i]
        frac = if abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) > (h/2)
            frac1 = (sign(dudt)*h/2Δt)/dudt
            frac2 = (sign(dvdt)*h/2Δt)/dvdt
            min(frac1, frac2)
        elseif abs(Δt*dudt) > (h/2) && abs(Δt*dvdt) < (h/2)
            (sign(dudt)*h/2Δt)/dudt
        elseif abs(Δt*dudt) < (h/2) && abs(Δt*dvdt) > (h/2)
            (sign(dvdt)*h/2Δt)/dvdt
        else
            1
        end
        if frac != 1
            @warn "Adjusting u and v velocities to prevent too high"
        end
        dudt = frac*dudt
        dvdt = frac*dvdt
        floes.u[i] += 1.5Δt*dudt-0.5Δt*floes.p_dudt[i]
        floes.v[i] += 1.5Δt*dvdt-0.5Δt*floes.p_dvdt[i]
        floes.p_dudt[i] = dudt
        floes.p_dvdt[i] = dvdt

        dξdt = (floes.trqOA[i] + ctrq)/floes.moment[i]
        dξdt = frac*dξdt
        ξ = floes.ξ[i] + 1.5Δt*dξdt-0.5Δt*floes.p_dξdt[i]
        if abs(ξ) > 1e-5
            @warn "Shrinking ξ" tstep = tstep
            ξ = sign(ξ) * 1e-5
        end
        floes.ξ[i] = ξ
        floes.p_dξdt[i] = dξdt

        # Update strain
        calc_strain!(LazyRow(floes, i))
    end
    return
end