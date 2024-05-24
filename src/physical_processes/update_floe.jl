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
        floe_settings,
        rng,
    )
Updates existing floe shape and related physical properties based of the polygon
defining the floe.
Inputs:
    floe        <Union{Floe, LazyRow{Floe}}> floe to update
    new_poly    <LG.Polygon> polygon representing new outline of floe
    new_mass    <AbstractFloat> mass of floe
    floe_settings   <FloeSettings> simulation's settings for making floes
    rng         <RNG> random number generator
Ouputs:
    Updates a given floe's physical properties given new shape and total mass.
"""
function replace_floe!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    new_poly,
    new_mass,
    floe_settings,
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
    floe.height = new_mass/(floe.area * floe_settings.ρi)
    floe.mass = new_mass
    floe.moment = calc_moment_inertia(
        floe.coords,
        floe.centroid,
        floe.height;
        ρi = floe_settings.ρi,
    )
    floe.angles = calc_poly_angles(floe.coords)
    floe.α = FT(0)
    translate!(floe.coords, -floe.centroid[1], -floe.centroid[2])
    floe.rmax = sqrt(maximum([sum(c.^2) for c in floe.coords[1]]))
    # Floe monte carlo points
    x_subfloe_points, y_subfloe_points, status = generate_subfloe_points(
        floe_settings.subfloe_point_generator,
        floe.coords,
        floe.rmax,
        floe.area,
        floe.status,
        rng,
    )
    translate!(floe.coords, floe.centroid[1], floe.centroid[2])
    floe.x_subfloe_points = x_subfloe_points
    floe.y_subfloe_points = y_subfloe_points
    # Floe status / identification
    floe.status = status
    return
end

"""
    conserve_momentum_change_floe_shape!(
        mass_tmp,
        moment_tmp,
        x_tmp,
        y_tmp,
        Δt,
        keep_floe,
        combine_floe = nothing,
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
    combine_floe <Union{Floe, LazyRow{Floe}}> if keep_floe's shape has been
                    changed due to an interaction with another floe, combine_floe
                    is that floe - optional parameter
Output:
    None. keep_floe's u, v, ξ, p_dxdt, p_dydt, p_dαdt, p_dudt, p_dvdt, and
    p_dξdt fields all updated to preserve momentum. 
Note:
    Function does not depend on conservation of mass
"""
function conserve_momentum_change_floe_shape!(
    mass_tmp,
    moment_tmp,
    x_tmp,
    y_tmp,
    Δt,
    keep_floe,
    combine_floe = nothing,
)
    # Calculate linear velocities to conserve linear momentum
    new_u = keep_floe.u * mass_tmp
    new_v = keep_floe.v * mass_tmp
    new_dxdt = keep_floe.p_dxdt * mass_tmp
    new_dydt = keep_floe.p_dydt * mass_tmp
    if !isnothing(combine_floe)
        new_u += combine_floe.u * combine_floe.mass
        new_v += combine_floe.v * combine_floe.mass
        new_dxdt += combine_floe.p_dxdt * combine_floe.mass
        new_dydt += combine_floe.p_dydt * combine_floe.mass
    end
    new_u /= keep_floe.mass
    new_v /= keep_floe.mass
    new_dxdt /= keep_floe.mass
    new_dydt /= keep_floe.mass
    # Calculate angular velocities to conserve rotational angular momentum
    p_x = x_tmp - Δt * keep_floe.p_dxdt
    p_y = y_tmp - Δt * keep_floe.p_dydt
    new_ξ = keep_floe.ξ * moment_tmp +  # spin momentum + orbital momentum
        mass_tmp * (x_tmp * keep_floe.v - y_tmp * keep_floe.u)
    new_dαdt = keep_floe.p_dαdt * moment_tmp +
        mass_tmp * (p_x * keep_floe.p_dydt - p_y * keep_floe.p_dxdt)
    if !isnothing(combine_floe)
        p_x = combine_floe.centroid[1] - Δt * combine_floe.p_dxdt
        p_y = combine_floe.centroid[2] - Δt * combine_floe.p_dydt
        new_ξ += combine_floe.ξ * combine_floe.moment +
            combine_floe.mass * (
                combine_floe.centroid[1] * combine_floe.v -
                combine_floe.centroid[2] * combine_floe.u
            )
        new_dαdt += combine_floe.p_dαdt * combine_floe.moment +
            combine_floe.mass * (
                p_x * combine_floe.p_dydt -
                p_y * combine_floe.p_dxdt
            )
    end
    # Subtract new orbital velocity from total previous momentum
    p_x = keep_floe.centroid[1] - Δt * new_dxdt
    p_y = keep_floe.centroid[2] - Δt * new_dydt
    new_ξ -= keep_floe.mass * (
        keep_floe.centroid[1] * new_v -
        keep_floe.centroid[2] * new_u
    )
    new_dαdt -= keep_floe.mass * (p_x * new_dydt - p_y * new_dxdt)
    new_ξ /= keep_floe.moment
    new_dαdt /= keep_floe.moment

    # Set new values
    keep_floe.u = new_u
    keep_floe.v = new_v
    keep_floe.ξ = new_ξ
    keep_floe.p_dxdt = new_dxdt
    keep_floe.p_dydt = new_dydt
    keep_floe.p_dαdt = new_dαdt
    # Calculate previous accelerations
    keep_floe.p_dudt = (keep_floe.u - keep_floe.p_dxdt) / Δt
    keep_floe.p_dvdt = (keep_floe.v - keep_floe.p_dydt) / Δt
    keep_floe.p_dξdt = (keep_floe.ξ - keep_floe.p_dαdt) / Δt
    return
end

"""
    update_new_rotation_conserve!(
        x, y,
        floe1, floe2,
        init_rot_momentum, init_p_rot_momentum,
        diff_orbital, diff_p_orbital,
        Δt,
    )

Update rotational velocities of two floes whose shapes may have changed to
maintain conservation of momentum given the floes' intial roational momentums
and the change in orbital momentum between their original shape and new shapes.
Additionally, both floes will have the same velocity at a given (x,y) point.
Inputs:
    x                   <Float> x-coordinate of point where floes share same
                            rotational velocity
    y                   <Float> y-coordinate of point where floes share same
                            rotational velocity
    floe1               <Union{Floe, LazyRow{Floe}}> first floe in pair
    floe2               <Union{Floe, LazyRow{Floe}}> second floe in pair
    init_rot_momentum   <Float> initial rotational momentum of floe 1 before
                            shape change
    init_p_rot_momentum <Float> initial rotational momentum of floe 2 before
                            shape change
    diff_orbital        <Float> change in floe 1's orbital velocity after shape
                            change
    diff_p_orbital      <Float> change in floe 1's orbital velocity after shape
                            change
    Δt                  <Int> length of timestep in seconds
Output:
    Nothing. Update both floes' rotational velocity, previous rotational
    velocity, and rotational acceleration to conserve angular momentum.
"""
function update_new_rotation_conserve!(
    x, y,
    floe1, floe2,
    init_rot_momentum, init_p_rot_momentum,
    diff_orbital, diff_p_orbital,
    Δt,
)
    # Find radius of each polygon to shared midpoint
    #x, y = find_shared_edges_midpoint(floe1.coords, floe2.coords)
    rad1 = sqrt(
        (floe1.centroid[1] - x)^2 +
        (floe1.centroid[2] - y)^2
    )
    rad2 = sqrt(
        (floe2.centroid[1] - x)^2 +
        (floe2.centroid[2] - y)^2
    )
    rad_ratio = rad1 / rad2
    # Determine ξ values so they are stationary at midpoint of shared edges
    floe1.ξ = (diff_orbital + init_rot_momentum) /
        (floe1.moment - (floe2.moment * rad_ratio))
    floe2.ξ = -floe1.ξ * rad_ratio
    # Determine p_dαdt values so they are stationary at intersection point
    floe1.p_dαdt = (diff_p_orbital + init_p_rot_momentum) /
        (floe1.moment - (floe2.moment * rad_ratio))
    floe2.p_dαdt = -floe1.p_dαdt * rad_ratio
    # Calculate previous rotational accelerations
    floe1.p_dξdt = (floe1.ξ - floe1.p_dαdt) / Δt
    floe2.p_dξdt = (floe2.ξ - floe2.p_dαdt) / Δt
end

"""
    conserve_momentum_fracture_floe!(
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
Note: Depends on conservation of mass.
"""
function conserve_momentum_fracture_floe!(
    init_floe,
    new_floes::StructArray{<:Floe{FT}},
    Δt,
) where {FT}
    if !isempty(new_floes)
        # conserve linear momentum by keeping linear velocities the same
        new_floes.u .= init_floe.u 
        new_floes.v .= init_floe.v
        new_floes.p_dxdt .= init_floe.p_dxdt
        new_floes.p_dydt .= init_floe.p_dydt
        new_floes.p_dudt .= init_floe.p_dudt
        new_floes.p_dvdt .= init_floe.p_dvdt
        # Reset α with new coordinates
        new_floes.α .= 0
        # MATLAB uses assumptions -> doesn't conserve rotational momentum
        new_floes.ξ .= init_floe.ξ
        new_floes.p_dαdt .= 0
        new_floes.p_dξdt .= init_floe.p_dξdt
    end
end

"""
    conserve_momentum_transfer_mass!(
        floes,
        idx1, idx2,
        m1, m2,
        I1, I2,
        x1, x2,
        y1, y2,
        Δt,
        pieces_list = nothing,
        pieces_idx = 0,
    )

Conserve linear momentum when mass is transfered from one floe to another floe.
Inputs:
    floes       <StructArray{Floes}> list of floes
    idx1        <Int> index of first floe in floes list
    idx2        <Int> index of second floe in floes list
    m1          <Float> initial mass of first floe before mass transfer
    m2          <Float> initial mass of second floe before mass transfer
    Δt          <Int> timestep of simulation in seconds
    pieces_list <Union{nothing, Vector{Floe}}> list of floes created from floes
                    breaking during ridging and rafting
    pieces_idx  <Int> start index of pieces included in conservation
                    calculations within the pieces_list
Outputs:
    Nothing. Update floes' velocities and accelerations to conserve linear
    momentum.
"""
function conserve_momentum_transfer_mass!(
    floes,
    idx1, idx2,
    m1, m2,
    Δt,
    pieces_list = nothing,
    pieces_idx = 0,
)
    # Conserve linear - assume resultant floes have same linear velocity
    tot_mass = m1 + m2
    new_u = (m1 * floes.u[idx1] + m2 * floes.u[idx2]) / tot_mass
    new_v = (m1 * floes.v[idx1] + m2 * floes.v[idx2]) / tot_mass
    new_p_dxdt = (m1 * floes.p_dxdt[idx1] + m2 * floes.p_dxdt[idx2]) / tot_mass
    new_p_dydt = (m1 * floes.p_dydt[idx1] + m2 * floes.p_dydt[idx2]) / tot_mass
    new_p_dudt = (floes.u[idx1] - floes.p_dxdt[idx1]) / Δt
    new_p_dvdt = (floes.v[idx1] - floes.p_dydt[idx1]) / Δt
    # MATLAB uses assumptions -> doesn't conserve rotational momentum
    floes.p_dαdt[idx1] = 0
    floes.p_dαdt[idx2] = 0
    # Update pieces list of they exist
    if !isnothing(pieces_list) && length(pieces_list) > 0
        pieces_list.p_dαdt[pieces_idx:end] .= 0
        pieces_list.u[pieces_idx:end] .= new_u
        pieces_list.v[pieces_idx:end] .= new_v
        pieces_list.p_dxdt[pieces_idx:end] .= new_p_dxdt
        pieces_list.p_dydt[pieces_idx:end] .= new_p_dydt
        pieces_list.p_dudt[pieces_idx:end] .= new_p_dudt
        pieces_list.p_dvdt[pieces_idx:end] .= new_p_dvdt
    end
    # Update floes linear velocities and accelerations
    floes.u[idx1], floes.u[idx2] = new_u, new_u
    floes.v[idx1], floes.v[idx2] = new_v, new_v
    floes.p_dxdt[idx1], floes.p_dxdt[idx2] = new_p_dxdt, new_p_dxdt
    floes.p_dydt[idx1], floes.p_dydt[idx2] = new_p_dydt, new_p_dydt
    floes.p_dudt[idx1], floes.p_dudt[idx2] = new_p_dudt, new_p_dudt
    floes.p_dvdt[idx1], floes.p_dvdt[idx2] = new_p_dvdt, new_p_dvdt
    return
end

"""
    update_ghost_timestep_vals!(floes, idx, parent_idx)

Update a parent floes and its ghosts to match velocities and accelerations at
given index.
Inputs:
    floes       <StructArray{Floe}> list of floes in the simulation
    idx         <Int> index of floe within floes list to copy to floe and ghosts
    parent_idx  <Int> index of parent floe to update, along with its ghosts
Outputs:
    Nothing. Update floe values. 
"""
function update_ghost_timestep_vals!(floes, idx, parent_idx)
    floes.u[parent_idx] = floes.u[idx]
    floes.v[parent_idx] = floes.v[idx]
    floes.ξ[parent_idx] = floes.ξ[idx]
    floes.p_dxdt[parent_idx] = floes.p_dxdt[idx]
    floes.p_dydt[parent_idx] = floes.p_dydt[idx]
    floes.p_dudt[parent_idx] = floes.p_dudt[idx]
    floes.p_dvdt[parent_idx] = floes.p_dvdt[idx]
    floes.p_dαdt[parent_idx] = floes.p_dαdt[idx]
    floes.p_dξdt[parent_idx] = floes.p_dξdt[idx]
    for gidx in floes.ghosts[parent_idx]
        floes.u[gidx] = floes.u[parent_idx]
        floes.v[gidx] = floes.v[parent_idx]
        floes.ξ[gidx] = floes.ξ[parent_idx]
        floes.p_dxdt[gidx] = floes.p_dxdt[parent_idx]
        floes.p_dydt[gidx] = floes.p_dydt[parent_idx]
        floes.p_dudt[gidx] = floes.p_dudt[parent_idx]
        floes.p_dvdt[gidx] = floes.p_dvdt[parent_idx]
        floes.p_dαdt[gidx] = floes.p_dαdt[parent_idx]
        floes.p_dξdt[gidx] = floes.p_dξdt[parent_idx]
    end
    return
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
function calc_stress!(floe::Union{LazyRow{Floe{FT}}, Floe{FT}}, floe_settings, Δt) where {FT}
    # Stress calcultions
    xi, yi = floe.centroid
    inters = floe.interactions
    # Calculates timestep stress
    stress = fill(FT(0), 2, 2)
    for i in 1:floe.num_inters
        stress[1, 1] += (inters[i, xpoint] - xi) * inters[i, xforce]
        stress[1, 2] += (inters[i, ypoint] - yi) * inters[i, xforce] +
            (inters[i, xpoint] - xi) * inters[i, yforce]
        stress[2, 2] += (inters[i, ypoint] - yi) * inters[i, yforce]
    end
    stress[1, 2] *= FT(0.5)
    stress[2, 1] = stress[1, 2]
    stress .*= 1/(floe.area * floe.height)
    update_damage!(floe_settings.stress_calculator, stress, floe, Δt)
    update_stress!(floe_settings.stress_calculator, stress, floe)
    # push!(floe.stress_instant, stress)
    floe.stress_instant = stress

    return
end

# Get rid of typing for decay calc. 
function update_damage!(stress_calculator, curr_stress, floe, Δt)
    return
end

# function update_damage!(different ytpe of stress calculator)
function update_damage!(stress_calculator::DamageStressCalculator, curr_stress, floe, Δt)
    τ = stress_calculator.τ
    if iszero(curr_stress)
        floe.damage = 0
    else
        floe.damage = (floe.stress_accum + (Δt/τ)*(curr_stress - floe.stress_accum))./curr_stress
    end
end

function update_stress!(stress_calculator::DamageStressCalculator, curr_stress, floe)
    floe.stress_accum = floe.damage * curr_stress
end

function update_stress!(stress_calculator, curr_stress, floe)
    floe.stress_accum = floe.stress_accum + floe.damage*(curr_stress - floe.stress_accum)
end

"""
    calc_strain!(coords, centroid, u, v, ξ, area)

Calculates the strain on a floe given the velocity at each vertex
Inputs:
    floe        <Floe{AbstractFloat}> a floe
Outputs:
    strain      <Matrix{AbstractFloat}> 2x2 matrix for floe strain 
"""
function calc_strain!(floe::Union{LazyRow{Floe{FT}}, Floe{FT}}) where {FT}
    # coordinates of floe centered at centroid
    translate!(floe.coords, -floe.centroid[1], -floe.centroid[2])
    fill!(floe.strain, FT(0))
    for i in 1:(length(floe.coords[1]) - 1)
        xdiff = floe.coords[1][i + 1][1] - floe.coords[1][i][1]
        ydiff = floe.coords[1][i + 1][2] - floe.coords[1][i][2]
        rad1 = sqrt(floe.coords[1][i][1]^2 + floe.coords[1][i][2]^2)
        θ1 = atan(floe.coords[1][i][2], floe.coords[1][i][1])
        rad2 = sqrt(floe.coords[1][i + 1][1]^2 + floe.coords[1][i + 1][2]^2)
        θ2 = atan(floe.coords[1][i + 1][2], floe.coords[1][i + 1][1])
        u1 = floe.u - floe.ξ * rad1 * sin(θ1)
        u2 = floe.u - floe.ξ * rad2 * sin(θ2)
        v1 = floe.u + floe.ξ * rad1 * cos(θ1)
        v2 = floe.u + floe.ξ * rad2 * cos(θ2)
        udiff = u2 - u1
        vdiff = v2 - v1
        floe.strain[1, 1] += udiff * ydiff
        floe.strain[1, 2] += udiff * xdiff + vdiff * ydiff
        floe.strain[2, 2] += vdiff * xdiff
    end
    floe.strain[1, 2] *= FT(0.5)
    floe.strain[2, 1] = floe.strain[1, 2]
    floe.strain ./= 2floe.area
    translate!(floe.coords, floe.centroid[1], floe.centroid[2])
    return
end

"""
    timestep_floe(floe)

Update floe position and velocities using second-order time stepping with
tendencies calculated at previous timesteps. Height, mass, stress, and strain
also updated based on previous timestep thermodynamics and interactions with
other floes. 
Input:
        floe            <Floe>
        Δt              <Int> simulation timestep in second
        floe_settings   <FloeSettings> simulation floe settings
Output:
        None. Floe's fields are updated with values.
"""
function timestep_floe_properties!(
    floes,
    tstep,
    Δt,
    floe_settings,
)
    Threads.@threads for i in eachindex(floes)
        cforce = floes.collision_force[i]
        ctrq = floes.collision_trq[i]
        # Update stress
        if floes.num_inters[i] > 0
            calc_stress!(LazyRow(floes, i), floe_settings, Δt)
        end
        # Ensure no extreem values due to model instability
        if floes.height[i] > floe_settings.max_floe_height
            @warn "Reducing height to 10 m"
            floes.height[i] = floe_settings.max_floe_height
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