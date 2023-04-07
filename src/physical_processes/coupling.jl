"""
Functions needed for coupling the ice, ocean, and atmosphere.
"""

#-------------- Monte Carlo Point Calculations --------------#

"""
    find_centered_cell_indices(xp, yp, grid::RegRectilinearGrid)

Find indicies of cells centered on grid lines of the given RegRectilinearGrid
that the given x-coordinates and y-coordinates fall within.
These cells are centered around the grid lines, so they are shifted grid cells
by half a cell. Method depends on grid being a regular rectilinear grid.
Inputs:
    xp      <Vector{Float}> x-coordinates of a vector of points
    yp      <Vector{Float}> y-coordinates of a vector of points
    grid    <RegRectilinearGrid> simulation grid
Outputs:
    xidx    <Vector{Float}> x-indicies of grid cell (cented on grid lines) each
                x point is within
    yidx    <Vector{Float}> y-indicies of grid cell (cented on grid lines) each
                y point is within
Note:
    Points can be outside of the grid, so indices can be less than 1 or greater
    than the number of grid lines in a given direction.
"""
function find_centered_cell_indices(xp, yp, grid::RegRectilinearGrid)
    xidx = floor.(Int, (xp .- grid.xg[1])/(grid.xg[2] - grid.xg[1]) .+ 0.5) .+ 1
    yidx = floor.(Int, (yp .- grid.yg[1])/(grid.yg[2] - grid.yg[1]) .+ 0.5) .+ 1
    return xidx, yidx
end

"""
    filter_oob_points(
        p,
        xr,
        yr,
        grid,
        ::NonPeriodicBoundary,
        ::NonPeriodicBoundary,
    )

With all non-periodic boundaries, points outside of the grid in both the x and y
direction need to be removed. These points can't be interpolated as we don't
have any information on the ocean outside of the grid. 
Note that p is a matrix of non-translated points centered on (0,0),
while xr and yr are these coordinates translated by a floe's centroid. 
Inputs:
    p   <Matrix{Float}> a 2xn matrix of points where the 1st row corresponds
            to x-values and the 2nd row corresponds to y-values
    xr  <Vector{Float}> a length-n vector of translated x-coordinates
    yr  <Vector{Float}> a length-n vector of translated y-coordinates
        <::NonPeriodicBoundary> type of either north or south boundary -
            checking if periodic pair
        <::NonPeriodicBoundary> type of either east or west boundary -
            checking if periodic pair
Output:
    p, xr, and yr filtered so that the translated points are all within the grid
"""
function filter_oob_points(
    p,
    xr,
    yr,
    grid,
    ::NonPeriodicBoundary,
    ::NonPeriodicBoundary,
)
    keep_idx = filter(i -> (grid.xg[1] .<= xr[i] .<= grid.xg[end]) && 
        (grid.yg[1] .<= yr[i] .<= grid.yg[end]), 1:length(xr))
    return p[:, keep_idx], xr[keep_idx], yr[keep_idx]
end

"""
    filter_oob_points(
        p,
        xr,
        yr,
        grid,
        ::NonPeriodicBoundary,
        ::PeriodicBoundary,
    )

With non-periodic boundaries in the north/south direction, points outside of the
grid in the y direction need to be removed. These points can't be interpolated
as we don't have any information on the ocean outside of the grid. In the
periodic direction, we can't have a floe longer than the grid, as the floe could
overlap with itself.
Note that p is a matrix of non-translated points centered on (0,0),
while xr and yr are these coordinates translated to a floe location. 
Inputs:
    p   <Matrix{Float}> a 2xn matrix of points where the 1st row corresponds to
            x-values and the 2nd row corresponds to y-values
    xr  <Vector{Float}> a length-n vector of translated x-coordinates
    yr  <Vector{Float}> a length-n vector of translated y-coordinates
        <::NonPeriodicBoundary> type of either north or south boundary -
            checking if periodic pair
        <::PeriodicBoundary> type of either east or west boundary - checking if
            periodic pair
Output:
    p, xr, and yr filtered so that the translated points are all within the grid
    in the y-direction and so floe does not overlap with itself through periodic
    boundary. 
"""
function filter_oob_points(
    p,
    xr,
    yr,
    grid,
    ::NonPeriodicBoundary,
    ::PeriodicBoundary,
)
    ymin, ymax = extrema(yr)
    keep_idx =
        if (ymax - ymin) > (grid.yg[end] - grid.yg[1])
            @warn "A floe longer than the domain passed through a periodic \
                boundary. It was removed to prevent overlap."
            Vector{Int}(undef, 0)
        else
            findall(p -> grid.xg[1] .<= p .<= grid.xg[end], xr)
        end
    return p[:, keep_idx], xr[keep_idx], yr[keep_idx]
end

"""
    filter_oob_points(
        p,
        xr,
        yr,
        grid,
        ::PeriodicBoundary,
        ::NonPeriodicBoundary,
    )

With non-periodic boundaries in the east/west direction, points outside of the
grid in the x direction need to be removed. These points can't be interpolated
as we don't have any information on the ocean outside of the grid. In the
periodic direction, we can't have a floe longer than the grid as the floe could
overlap with itself.
Note that p is a matrix of non-translated points centered on (0,0), while xr and
yr are these coordinates translated to a floe location. 
Inputs:
    p   <Matrix{Float}> a 2xn matrix of points where the 1st row corresponds to
            x-values and the 2nd row corresponds to y-values
    xr  <Vector{Float}> a length-n vector of translated x-coordinates
    yr  <Vector{Float}> a length-n vector of translated y-coordinates
        <::PeriodicBoundary> type of either north or south boundary - checking
            if periodic pair
        <::AbstractBoundary> type of either east or west boundary - checking if
            periodic pair
Output:
    p, xr, and yr filtered so that the translated points are all within the grid
    in the x-direction and so that the floe does not overlap with itself through
    periodic boundary. 
"""
function filter_oob_points(
    p,
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::NonPeriodicBoundary,
)
    xmin, xmax = extrema(xr)
    keep_idx =
        if (xmax - xmin) > (grid.xg[end] - grid.xg[1])
            @warn "A floe longer than the domain passed through a periodic \
                boundary. It was removed to prevent overlap."
            Vector{Int}(undef, 0)
        else
            findall(p -> grid.yg[1] .<= p .<= grid.yg[end], yr)
        end
    return p[:, keep_idx], xr[keep_idx], yr[keep_idx]
end

"""
    filter_oob_points(
        p,
        xr,
        yr,
        grid,
        ::PeriodicBoundary,
        ::PeriodicBoundary,
    )

With all periodic boundaries, no points should be filtered, unless floe is
longer than the grid as the floe could overlap with itself. In this case, we
return no points. p is a matrix of non-translated points centered on (0,0),
while xr and yr are these coordinates translated to a floe location. 
Inputs:
    p   <Matrix{Float}> a 2xn matrix of points where the 1st row corresponds to
            x-values and the 2nd row corresponds to y-values
    xr  <Vector{Float}> a length-n vector of translated x-coordinates
    yr  <Vector{Float}> a length-n vector of translated y-coordinates
        <::PeriodicBoundary> type of either north or south boundary - checking
            if periodic pair
        <::PeriodicBoundary> type of either east or west boundary - checking if
            periodic pair
Output:
    p, xr, and yr as given, unless the floe might overlap with itself through
    periodic boundaries, in which case no points are returned.
"""
function filter_oob_points(
    p,
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::PeriodicBoundary,
)
    xmin, xmax = extrema(xr)
    ymin, ymax = extrema(yr)
    keep_idx =
        if (xmax - xmin) > (grid.xg[end] - grid.xg[1]) ||
            (ymax - ymin) > (grid.yg[end] - grid.yg[1])
                @warn "A floe longer than the domain passed through a periodic \
                    boundary. It was removed to prevent overlap."
                Vector{Int}(undef, 0)
        else
            collect(1:length(xr))
        end
    return p[:, keep_idx], xr[keep_idx], yr[keep_idx]
end

"""
    calc_mc_values(
        floe,
        grid,
        domain,
        ::Type{T} = Float64,
    )

Calculates monte carlo point's cartesian coordiantes, polar coordiantes,
velocity and index within the grid. 
Inputs:
    floe    <Floe> single floe
    grid    <AbstractGrid> model grid
    domain  <Domain> model domain
    T       <Type(AbstractFloat)> type to run simulation calculations
Outputs:
    mc_cart     <Matrix{AbstractFloat}> cartesian coordinates for model
                    coordinates - nx2 matrix of monte carlo coordinates where
                    first column is the x-coords and the second column is the
                    y-coords
    mc_polar    <Matrix{AbstractFloat}> polar coordinates for model coordinates
                    - nx2 matrix of monte carlo coordinates where first column
                    is the radius and the second column is the angle
    mc_vels     <Matrix{AbstractFloat}> velocity at each monte carlo point - nx2
                    matrix where the first column is the u-velocity and the
                    second is the v-velocity for n monte carlo points
    mc_grid_idx <Matrix{Int}> index of monte carlo points within the grid - nx2
                    matrix of indices where the first column is the grid column
                    index and the second column is the grid row index
"""
function calc_mc_values(floe, grid, domain, ::Type{T} = Float64) where {T}
    # Translate/rotate monte carlo points to floe location/orientation
    α = floe.α
    α_rot = [cos(α) -sin(α); sin(α) cos(α)]
    mc_p = α_rot * [floe.mc_x'; floe.mc_y']
    mc_xr = mc_p[1, :] .+ floe.centroid[1]
    mc_yr = mc_p[2, :] .+ floe.centroid[2]

    # Filter points outside of non-periodic boundaries and check if floe
    # overlaps through periodic bounds
    mc_p, mc_xr, mc_yr = filter_oob_points(
        mc_p,
        mc_xr,
        mc_yr,
        grid,
        domain.east,
        domain.north,
    )
    mc_polar = zeros(T, length(mc_xr), 2)
    mc_vels = zeros(T, length(mc_xr), 2)
    mc_grid_idx = zeros(Int, length(mc_xr), 2)
    if !isempty(mc_p)
        # Monte carlo points polar coordinates
        mc_rad = sqrt.(mc_p[1, :].^2 .+ mc_p[2, :].^2)
        mc_θ = atan.(mc_p[2, :], mc_p[1, :])
        mc_polar[:, 1] .= mc_rad
        mc_polar[:, 2] .= mc_θ
        # Monte carlo point velocities
        mc_u = floe.u .- floe.ξ * mc_rad .* sin.(mc_θ)
        mc_v = floe.v .+ floe.ξ * mc_rad .* cos.(mc_θ)
        mc_vels[:, 1] .= mc_u
        mc_vels[:, 2] .= mc_v
        # Find point indices for cells centered on grid lines
        mc_cols, mc_rows = find_cell_indices(mc_xr, mc_yr, grid)
        mc_grid_idx[:, 1] .= mc_cols
        mc_grid_idx[:, 2] .= mc_rows
    end
    return hcat(mc_xr, mc_yr),  # cartesian coords
        mc_polar,  # polar coords
        mc_vels,  # velocities
        mc_grid_idx  # grid cell indices
end

#-------------- Interpolation of Ocean and Atmosphere --------------#

"""
    find_interp_knots(
        point_idx,
        glines,
        Δd::Int,
        ::PeriodicBoundary,
        ::Type{T} = Float64,
    )

Find indicies in list of grid lines that surround points with indicies
'point_idx', with a buffer of Δd indices on each side of the points. In this
case, the points are being considered near a periodic boundary, which means that
they can loop around to the other side of the grid. If these points exist, we
extend the grid lines to cover the points and buffer. 
Inputs:
    point_idx   <Vector{Int}> vector of point indices representing the grid line
                    they are nearest
    glines      <Vector{Float}> vector of all grid line values 
    Δd          <Int> number of buffer grid cells to include on either side of
                    the provided indicies 
                <PeriodicBoundary> dispatching on periodic boundary
                <Float> datatype to run simulation - either Float32 or Float64
Outputs:
    knots       <Vector{AbstractFloat}> interpolation knots - grid line values
    knot_idx    <Vector{Int}> - indices of grid line values within grid list of
                    grid lines.
Note:
    The grid values are extended if points expand past gridlines, however, the
    indices are within the grid. For example, consider a grid where the maximum
    grid value is 1e5, with grid cells of length 1e4. One of the knot values
    might be 1.1e5, however, its index would be 2 since the grid line at 1e5 is
    equivalent to the first grid line since it is periodic, and 1.1e5 is one
    grid cell length past that value.

    This function depends on the ocean being periodic in the given direction.
    We assume that first grid line and the last grid line are the same, and have
    the same values within the ocean/atmosphere. These are not repeated in the
    knots, but rather only one is used. So if there are 10 grid lines, grid line
    1 and 10 are the equivalent and we use grid line 1 exclusively. 
"""
function find_interp_knots(
    point_idx,
    glines,
    Δd::Int,
    ::PeriodicBoundary,
    ::Type{T} = Float64,
) where T 
    min_line, max_line = extrema(point_idx)
    nlines = length(glines)
    ncells = nlines - 1

    # Grid lines surrounding points with buffers
    # Point close to ith grid line could be between the i and i-1 grid line
    min_line -= (Δd + 1)
    # Point close to ith grid line could be between the i and i+1 grid line
    max_line += (Δd + 1)

    # Find out-of-bounds (oob) indices and the in-bounds (within grid) indices
     # Out of bounds on south or west side of domain 
    low_oob_idx = Vector{Int}(undef, 0)
    # Out of bounds on north or east side of domain 
    high_oob_idx = Vector{Int}(undef, 0)
    in_bounds_idx = 
        # Last gird line is equal to first grid line
        if min_line < 1 && max_line > ncells
            low_oob_idx = (min_line + ncells):ncells
            high_oob_idx = 1:(max_line-ncells)
            1:ncells
        elseif min_line < 1
            low_oob_idx = (min_line + ncells):ncells
            1:max_line
        elseif max_line > ncells
            high_oob_idx = 1:(max_line-ncells)
            min_line:ncells
        else
            min_line:max_line
        end
    # Combine above indices for knot indicies
    knot_idx = [low_oob_idx; in_bounds_idx; high_oob_idx]
    #= Adjust out-of-bound values by grid length so there isn't a jump in
    interpolation spacing =#
    Δg = glines[end] - glines[1]
    knots = [glines[low_oob_idx] .- Δg;
        glines[in_bounds_idx]; glines[high_oob_idx] .+ Δg]
    return knots, knot_idx
end

"""
    find_interp_knots(
        point_idx,
        glines,
        Δd::Int,
        ::NonPeriodicBoundary,
        ::Type{T} = Float64,
        )

Find indicies in list of grid lines that surround points with indicies
'point_idx' with a buffer of Δd indices on each side of the points. In this
case, the points are being considered near a NON-periodic boundary, so we cut
off the possible indices past the edge of the grid. 
Inputs:
    point_idx   <Vector{Int}> vector of indices representing the grid line they
                    are nearest
    glines      <Vector{Float}> vector of grid line values 
    Δd          <Int> number of buffer grid cells to include on either side of
                    the provided indicies 
                <PeriodicBoundary> dispatching on periodic boundary
                <Float> datatype to run simulation - either Float32 or Float64
Outputs:
    knots       <Vector{AbstractFloat}> interpolation knots - grid line values
    knot_idx    <Vector{Int}> - indices of grid line values within grid list of
                    grid lines.
Note:
    Only knots within the grid will be returned since this is a non-periodic
    boundary.
"""
function find_interp_knots(
    point_idx,
    glines,
    Δd::Int,
    ::NonPeriodicBoundary,
    ::Type{T} = Float64,
) where T
    nlines = length(glines)
    min_line, max_line = extrema(point_idx)
    # point close to ith grid line could be between the i and i-1 grid line
    min_line -= (Δd + 1)
    # point close to ith grid line could be between the i and i+1 grid line
    max_line += (Δd + 1)
    # Can't be outside grid given unknown ocean/atmosphere conditions
    min_line = (min_line < 1) ? 1 : min_line
    max_line = (max_line > nlines) ? nlines : max_line
    return glines[min_line:max_line], [min_line:max_line;]
end

"""
    mc_interpolation(
        mc_cart,
        mc_grid_idx,
        grid,
        domain,
        atmos,
        ocean,
        coupling_settings,
        ::Type{T} = Float64
    )

Interpolates the ocean, atmosphere, and heatflux factor fields onto a floe's
monte carlo points.
Inputs:
    mc_cart             <Matrix{AbstractFloat}> cartesian coordinates for model
                            coordinates - nx2 matrix of monte carlo coordinates
                            where first column is the x-coords and the second
                            column is the y-coords
    mc_grid_idx         <Matrix{Int}> index of monte carlo points within the
                            grid - nx2 matrix of indices where the first column
                            is the grid column index and the second column is
                            the grid row index
    grid                <AbstractGrid> model grid
    domain              <Domain> model domain
    atmos               <Atmos> model atmosphere
    ocean               <Ocean> model ocean
    coupling_settings   <CouplingSettings> simulation's coupling settings
    T                   <Type(AbstractFloat)> type to run simulation calculations
Outputs:
    uatm        <Vector{AbstractFloat}> atmosphere u velocity interpolated onto
                    floe's monte carlo points
    vatm        <Vector{AbstractFloat}> atmosphere v velocity interpolated onto
                    floe's monte carlo points
    uocn        <Vector{AbstractFloat}> ocean u velocity interpolated onto
                    floe's monte carlo points
    vocn        <Vector{AbstractFloat}> ocean v velocity interpolated onto
                    floe's monte carlo points
    hflx_factor <Vector{AbstractFloat}> ocean heatflux factor interpolated onto
                    floe's monte carlo points
"""
function mc_interpolation(
    mc_cart,
    mc_grid_idx,
    grid,
    domain,
    atmos,
    ocean,
    coupling_settings,
    ::Type{T} = Float64
) where T<:AbstractFloat
    mc_xr = mc_cart[:, 1]
    mc_yr = mc_cart[:, 2]
    mc_cols = mc_grid_idx[:, 1]
    mc_rows = mc_grid_idx[:, 2]

    # Find knots and indices of knots for monte carlo interpolation
    xknots, xknot_idx = find_interp_knots(
        mc_cols,
        grid.xg,
        coupling_settings.Δd,
        domain.east,
        T,
    )
    yknots, yknot_idx = find_interp_knots(
        mc_rows,
        grid.yg,
        coupling_settings.Δd,
        domain.north,
        T,
    )
    # Atmos Interpolation for Monte Carlo Points
    uatm_interp = linear_interpolation(
        (yknots, xknots),
        atmos.u[yknot_idx, xknot_idx],
    )
    vatm_interp = linear_interpolation(
        (yknots, xknots),
        atmos.v[yknot_idx, xknot_idx],
    )
    uatm = [uatm_interp(mc_yr[i], mc_xr[i]) for i in eachindex(mc_xr)]
    vatm = [vatm_interp(mc_yr[i], mc_xr[i]) for i in eachindex(mc_xr)]

    # Ocean Interpolation for Monte Carlo Points
    uocn_interp = linear_interpolation(
        (yknots, xknots),
        ocean.u[yknot_idx, xknot_idx],
    )
    vocn_interp = linear_interpolation(
        (yknots, xknots),
        ocean.v[yknot_idx, xknot_idx],
    )
    uocn = [uocn_interp(mc_yr[i], mc_xr[i]) for i in eachindex(mc_xr)]
    vocn = [vocn_interp(mc_yr[i], mc_xr[i]) for i in eachindex(mc_xr)]

    # Heatflux Factor Interpolatoin for Monte Carlo Points
    hflx_interp = linear_interpolation(
        (yknots, xknots),
        ocean.hflx_factor[yknot_idx, xknot_idx],
    )
    hflx_factor = [hflx_interp(mc_yr[i], mc_xr[i]) for i in eachindex(mc_xr)]

    return uatm, vatm, uocn, vocn, hflx_factor
end

#-------------- Effects of Ice and Atmosphere on Ocean --------------#
"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::PeriodicBoundary,
        ::PeriodicBoundary,
    )

Return cell bounding values as is given the domain is doubley periodic and thus
the cell can extend beyond the grid as it will simply wrap back around into grid
through opposite periodic boundary.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractGrid> model's grid
            <PeriodicBoundary> type of north or south boundary - periodic pair
            <PeriodicBoundary> type of east or west boundary - periodic pair
Output:
    x and y minimums and maximums as given since they can extend past the grid
    due to periodic boundaries.
"""
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::PeriodicBoundary,
    ::PeriodicBoundary,
)
    return xmin, xmax, ymin, ymax
end

"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::NonPeriodicBoundary,
        ::PeriodicBoundary,
    )

Trim cell bound in the north-south direction if it exends past grid due to
non-periodic boundary pair.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractGrid>
            <NonPeriodicBoundary> type of either north or south boundary - not a
                periodic pair
            <PeriodicBoundary> type of either east or west boundary - periodic
                pair
Output:
        Potentially trimmed y min and y max if these values extend beyond grid
        values. Else returned unchanged. 
"""
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::NonPeriodicBoundary,
    ::PeriodicBoundary,
)
    ymin = ymin < grid.yg[1] ?
        grid.yg[1] :
        (ymin > grid.yg[end] ? grid.yg[end] : ymin)

    ymax = ymax > grid.yg[end] ?
        grid.yg[end] :
        (ymax < grid.yg[1] ? grid.yg[1] : ymax)
    return xmin, xmax, ymin, ymax
end

"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::PeriodicBoundary,
        ::NonPeriodicBoundary,
    )

Trim cell bound in the east-west direction if it exends past grid due to
non-periodic boundary pair.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractGrid>
            <PeriodicBoundary> type of either north or south boundary - periodic
                pair
            <NonPeriodicBoundary> type of either east or west boundary - not a
                periodic pair
Output:
    Potentially trimmed x min and max if these values extend beyond grid values.
    Else returned unchanged. 
"""
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::PeriodicBoundary,
    ::NonPeriodicBoundary,
)
    xmin = xmin < grid.xg[1] ?
        grid.xg[1] :
        (xmin > grid.xg[end] ? grid.xg[end] : xmin)

    xmax = xmax > grid.xg[end] ?
        grid.xg[end] :
        (xmax < grid.xg[1] ? grid.xg[1] : xmax) 
    return xmin, xmax, ymin, ymax
end

"""
    check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ::NonPeriodicBoundary,
        ::NonPeriodicBoundary,
    )

Trim cell bounds in the east-west and north-south direction if they exend past
grid due to non-periodic boundary pairs.
Inputs:
    xmin    <Float> center cell minimum x value
    xmax    <Float> center cell maxumum x value
    ymin    <Float> center cell minimum y value
    ymax    <Float> center cell maximum y value
    grid    <AbstractGrid>
            <NonPeriodicBoundary> type of either north or south boundary - not a
                periodic pair
            <NonPeriodicBoundary> type of either east or west boundary - not a
                periodic pair
Output:
    Potentially trimmed x and y minimums and maximums if these values extend
    beyond grid values. Else returned unchanged. 
"""
function check_cell_bounds(
    xmin,
    xmax,
    ymin,
    ymax,
    grid,
    ::NonPeriodicBoundary,
    ::NonPeriodicBoundary,
)
    xmin = xmin < grid.xg[1] ?
        grid.xg[1] :
        (xmin > grid.xg[end] ? grid.xg[end] : xmin)

    xmax = xmax > grid.xg[end] ?
        grid.xg[end] :
        (xmax < grid.xg[1] ? grid.xg[1] : xmax)

    ymin = ymin < grid.yg[1] ?
        grid.yg[1] :
        (ymin > grid.yg[end] ? grid.yg[end] : ymin)

    ymax = ymax > grid.yg[end] ?
        grid.yg[end] :
        (ymax < grid.yg[1] ? grid.yg[1] : ymax)
    return xmin, xmax, ymin, ymax
end

"""
    center_cell_coords(
        xidx::Int,
        yidx::Int,
        grid::RegRectilinearGrid,
        ns_bound,
        ew_bound,
    )

Find the coordinates of a given grid cell, centered on a grid line with row yidx
and column xidx. This is offset from the cells within the regular rectilinear
grid by half of a grid cell. 
Inputs:
    xidx        <Int> x index of grid line within list of gridlines (cell column)
    yidx        <Int> y index of grid line within list of gridlines (cell row)
    grid        <RegRectilinearGrid> model's grid 
    ns_bound    <AbstractBoundary> type of either north or south boundary - for
                    checking if periodic
    ew_bound    <AbstractBoundary> type of either east or west boundary - for
                    checking if perioidic
Output:
    <PolyVec> coordinates for cell centered on grid line with given indices.
    Note that cell bounds will be adjusted depending on if the bounds are
    periodic. Cells cannot extend outside of non-periodic boundaries and thus
    will be trimmed at boundaries. Therefore, if indices place cell completely
    outside of grid, could return a line at the edge of the boundary. 
"""
function center_cell_coords(
    xidx::Int,
    yidx::Int,
    grid::RegRectilinearGrid,
    ns_bound,
    ew_bound,
)
    Δx = grid.xg[2] .- grid.xg[1]
    Δy = grid.yg[2] .- grid.yg[1]
    xmin = (xidx - 1.5)*Δx + grid.xg[1]
    xmax = xmin + Δx
    ymin = (yidx - 1.5)*Δy + grid.yg[1]
    ymax = ymin + Δy
    #= Check if cell extends beyond boundaries and if non-periodic, trim cell to
    fit within grid. =#
    xmin, xmax, ymin, ymax = check_cell_bounds(
        xmin,
        xmax,
        ymin,
        ymax,
        grid,
        ns_bound,
        ew_bound,
    )
    return [[[xmin, ymin], [xmin, ymax],
    [xmax, ymax], [xmax, ymin],
    [xmin, ymin]]]
end

"""
    shift_cell_idx(idx, nlines, ::NonPeriodicBoundary)

Return index as is given non-periodic boundary pair in either x or y direction.
Inputs:
    idx     <Int> grid line index in either x or y
    nlines  <Int> number of grid lines in model grid in either x or y direction
            <NonPeriodicBoundary> boundary pair is non-periodic
Ouput:
    idx <Int> as given. Can include the index nlines, unlike with the periodic
        case, which will use the first index instead. 
"""
function shift_cell_idx(idx, nlines, ::NonPeriodicBoundary)
    return idx
end

"""
    shift_cell_idx(idx, nlines, ::PeriodicBoundary)

If index is greater than or equal to the grid lines, shift index to equivalent
grid line on opposite side of grid due to periodic boundary. Similarly if given
index is less than 1, shift index to equivalent grid line on opposite side of
grid due to periodic boundary.
Inputs:
    idx     <Int> grid line index in either x or y
    nlines  <Int> number of grid lines in model grid in either x or y direction
            <PeriodicBoundary> boundary pair is periodic
Output:
    <Int> if given index is greater than or equal to number of grid lines, shift
    index. If given index is less than 1, shift grid index. For example, the
    last grid index, nlines, is equivalent to the 1st grid line. The nlines+1
    grid line is equivalent to the 2nd grid line.
"""
function shift_cell_idx(idx, nlines, ::PeriodicBoundary)
    ncells = nlines - 1
    return idx < 1 ? (idx + ncells) : ncells < idx ? (idx - ncells) : idx
end

"""
    aggregate_grid_force!(
        mc_grid_idx,
        τx_ocn,
        τy_ocn,
        floe,
        ocean,
        grid,
        ns_bound,
        ew_bound,
        spinlock,
        ::Type{T} = Float64
    )

Add force from the ice on ocean to ocean force fields (fx & fy) for each grid
cell and update ocean sea ice area fraction (si_area), representing total area
of sea ice in a given cell.
Inputs:
    mc_grid_idx <Matrix{Int}> index of monte carlo points within the grid - nx2
                    matrix of indices where the first column is the grid column
                    index and the second column is the grid row index
    τx_ocn      <Vector{Float}> x-stress caused by each corresponding monte
                    carlo point on the ocean
    τy_ocn      <Vector{Float}> y-stress caused by each corresponding monte
                    carlo point on the ocean
    floe        <Floe> single floe within model
    ocean       <Ocean> model's ocean
    grid        <AbstractGrid> model's grid
    ns_bound    <AbstractBoundary> north or south boundary of domain -
                    dispatches on periodic vs non-periodic
    ew_bound    <AbstractBoundary> east or west boundary of domain -
                    dispatches on periodic vs non-periodic
    spinlock    <Thread.SpinLock>
    T           <Float> Float type simulation is using for calculations
                (Float32 or Float64)
Outputs:
    None. Ocean fields updated in-place. Note that the needed additions take
    place within a lock. Since multiple floes can be within one grid cell, this
    critical section needs to be locked for when the code is run with multiple
    threads to prevent race conditions.
"""
function aggregate_grid_force!(
    mc_grid_idx,
    τx_ocn,
    τy_ocn,
    floe,
    ocean,
    grid,
    ns_bound,
    ew_bound,
    spinlock, # lock for parallelization
    ::Type{T} = Float64
) where {T}
    mc_cols = mc_grid_idx[:, 1]
    mc_rows = mc_grid_idx[:, 2]
    τx_group = Dict{Tuple{Int64, Int64}, Vector{T}}()
    τy_group = Dict{Tuple{Int64, Int64}, Vector{T}}()
    # Use dictionary to sort stress values by grid cell
    for i in eachindex(mc_cols)
        idx = (mc_rows[i], mc_cols[i])
        if haskey(τx_group, idx)
            push!(τx_group[idx], τx_ocn[i])
            push!(τy_group[idx], τy_ocn[i])
        else
            τx_group[idx] = [τx_ocn[i]]
            τy_group[idx] = [τy_ocn[i]]
        end
    end

    # Determine force from floe on each grid cell it is in
    floe_poly = LG.Polygon(floe.coords)
    for (row, col) in keys(τx_group)
        grid_poly = LG.Polygon(center_cell_coords(
                col,
                row,
                grid,
                ns_bound,
                ew_bound
        ))
        floe_area_in_cell = LG.area(LG.intersection(grid_poly, floe_poly))
        if floe_area_in_cell > 0
            fx_mean = mean(τx_group[row, col]) * floe_area_in_cell
            fy_mean = mean(τy_group[row, col]) * floe_area_in_cell
            row = shift_cell_idx(row, grid.dims[1] + 1, ns_bound)
            col = shift_cell_idx(col, grid.dims[2] + 1, ew_bound)
            Threads.lock(spinlock) do
                ocean.fx[row, col] += fx_mean
                ocean.fy[row, col] += fy_mean
                ocean.si_area[row, col] += floe_area_in_cell
            end
        end
    end
    return
end

#-------------- Ocean and Atmosphere on Ice --------------#
"""
    calc_atmosphere_forcing(
        uatm,
        vatm,
        c,  # constants
    )

Calculates the stresses on a floe from the atmosphere above given the atmosphere
velocity at each of its monte carlo points.
Inputs:
    uatm    <Vector{AbstractFloat}> atmosphere u-velocity at each monte carlo
                point
    vatm    <Vector{AbstractFloat}> atmosphere v-velocity at each monte carlo
                point
    c       <Constants> simulation constants
Outputs:
    τx_atm  <Vector{AbstractFloat}> stress from atmosphere on floe in
                x-direction at each monte carlo point
    τy_atm  <Vector{AbstractFloat}> stress from atmosphere on floe in
                y-direction at each monte carlo point
"""
function calc_atmosphere_forcing(
    uatm,
    vatm,
    c,  # constants
)
    # Average velocity over the floe
    avg_uatm = mean(uatm)
    avg_vatm = mean(vatm)

    # Stress on ice from atmopshere
    τx_atm = (c.ρa * c.Cd_ia * sqrt(avg_uatm^2 + avg_vatm^2) * avg_uatm)
    τy_atm = (c.ρa * c.Cd_ia * sqrt(avg_uatm^2 + avg_vatm^2) * avg_vatm)
    return τx_atm, τy_atm
end

"""
    calc_ocean_forcing(
        floe,
        mc_vel,
        uocn,
        vocn,
        c,
    ) 

Calculates the stresses on a floe from the ocean below given the ocean velocity
at each of its monte carlo points.
Inputs:
    floe    <Floe> single model floe
    mc_vel  <Matrix{AbstractFloat}> velocity at each monte carlo point - nx2
                matrix where the first column is the u-velocity and the second
                is the v-velocity for n monte carlo points within floe
    uocn    <Vector{AbstractFloat}> ocean u-velocity at each monte carlo point
    vocn    <Vector{AbstractFloat}> ocean v-velocity at each monte carlo point
    c       <Constants> simulation constants
Outputs:
    τx_ocn          <Vector{AbstractFloat}> stress from ocean on floe in
                        x-direction at each monte carlo point
    τy_ocn          <Vector{AbstractFloat}> stress from ocean on floe in
                        y-direction at each monte carlo point
    τx_pressure∇    <Vector{AbstractFloat}> stress from ocean pressure gradient
                        in x-direction at each monte carlo point
    τy_pressure∇    <Vector{AbstractFloat}> stress from ocean pressure gradient
                        in y-direction at each monte carlo point
"""
function calc_ocean_forcing(
    floe,
    mc_vel,
    uocn,
    vocn,
    c,  # constants
)

    # Stress on ice from ocean
    Δu_OI = uocn .- mc_vel[:, 1]
    Δv_OI = vocn .- mc_vel[:, 2]
    τx_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .*
        (cos(c.turnθ) .* Δu_OI .- sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo*c.Cd_io*sqrt.(Δu_OI.^2 + Δv_OI.^2) .*
        (sin(c.turnθ) .* Δu_OI .+ cos(c.turnθ) * Δv_OI)

    # Stress on ice from pressure gradient
    ma_ratio = floe.mass/floe.area
    τx_pressure∇ = -ma_ratio * c.f .* vocn
    τy_pressure∇ = ma_ratio * c.f .* uocn

    return τx_ocn, τy_ocn, τx_pressure∇, τy_pressure∇
end

"""
    timestep_coupling!(
        model,
        consts,
        coupling_settings,
        spinlock,
        ::Type{T} = Float64
    )

Calculates the effects of the ocean and atmosphere on the ice and the effects of
the ice on the ocean.
Inputs:
    model               <Model> model
    consts              <Constants> constants used in simulation
    coupling_settings   <CouplingSettings> settings for coupling
    spinlock            <Thread.SpinLock>
    T                   <AbstractFloat> type to run simulation calculations
Outputs:
    None. Updates each floe's ocean/atmosphere forcings (fxOA, fyOA, torqueOA)
    and calculates stresses on each ocean grid cell.          
"""
function timestep_coupling!(
    model,
    consts,
    coupling_settings,
    spinlock,
    ::Type{T} = Float64,
) where T<:AbstractFloat
    # Clear ocean forces and area fractions
    fill!(model.ocean.fx, T(0))
    fill!(model.ocean.fy, T(0))
    fill!(model.ocean.si_area, T(0))
    # Calcualte coupling for each floe
    Threads.@threads for i in eachindex(model.floes)
        ifloe = model.floes[i]
        # Find monte carlo point peroperties
        mc_cart, mc_polar, mc_vel, mc_grid_idx = calc_mc_values(
            ifloe,
            model.grid,
            model.domain,
        )
        # Interpolate ocean and atmosphere values onto monte carlo points
        uatm, vatm, uocn, vocn, hflx_factor = mc_interpolation(
            mc_cart,
            mc_grid_idx,
            model.grid,
            model.domain,
            model.atmos,
            model.ocean,
            coupling_settings,
            T,
        )
        # Calculate effects on the floe
        ifloe.hflx_factor = mean(hflx_factor)
        τx_atm, τy_atm = calc_atmosphere_forcing(uatm, vatm, consts)
        τx_ocn, τy_ocn, τx_pressure∇, τy_pressure∇ = calc_ocean_forcing(
            ifloe,
            mc_vel,
            uocn,
            vocn,
            consts,
        )

        # Update ocean with force from floe per grid cell
        if coupling_settings.calc_ocnτ_on
            aggregate_grid_force!(
                mc_grid_idx,
                -τx_ocn,
                -τy_ocn,
                ifloe,
                model.ocean,
                model.grid,
                model.domain.east,
                model.domain.north,
                spinlock,
                T,
            )
        end

        # Sum stresses and find torques
        mc_rad = mc_polar[:, 1]
        mc_θ = mc_polar[:, 2]
        τx = τx_atm .+ τx_pressure∇ .+ τx_ocn
        τy = τy_atm .+ τy_pressure∇ .+ τy_ocn
        τtrq = (-τx .* sin.(mc_θ) .+ τy .* cos.(mc_θ)) .* mc_rad

        # Add coriolis force to total foces - does not contribute to torque
        ma_ratio = (ifloe.mass/ifloe.area)
        τx .+= ma_ratio * consts.f * ifloe.v
        τy .-= ma_ratio * consts.f * ifloe.u

        # Average forces on ice floe
        ifloe.fxOA = mean(τx) * ifloe.area
        ifloe.fyOA = mean(τy) * ifloe.area
        ifloe.trqOA = mean(τtrq) * ifloe.area
        model.floes[i] = ifloe
    end
end

"""
    timestep_ocean!(m, c, Δt)

Update model's ocean from affects of atmsophere and ice
Input:
        m   <Model> simulation model
        c   <Constants> simulation constants
        Δt  <Int> simulation timestep
Outputs: 
        None. The ocean stress fields are updated from atmos stress.
        Heatflux is also updated with ocean and atmos temperatures. 
"""
function timestep_ocean!(m, c, Δt)
    # Fraction of grid cells covered in ice vs ocean
    cell_area = (m.grid.xg[2] - m.grid.xg[1]) * (m.grid.yg[2] - m.grid.yg[1])
    si_frac = m.ocean.si_area ./ cell_area
    icy_idx = si_frac .> 0
    ocn_frac = 1 .- si_frac
    # Add ice stress on ocean
    m.ocean.τx .= m.ocean.fx
    m.ocean.τy .= m.ocean.fy
    m.ocean.τx[icy_idx] ./= m.ocean.si_area[icy_idx]
    m.ocean.τy[icy_idx] ./= m.ocean.si_area[icy_idx]
    # Atmospheric stress on ocean - resets ocean stresses
    Δu_AO = m.atmos.u .- m.ocean.u
    Δv_AO = m.atmos.v .- m.ocean.v
    vel_norm = sqrt.(Δu_AO.^2 + Δv_AO.^2)
    m.ocean.τx .+= c.ρa * c.Cd_ao * ocn_frac .* vel_norm .* Δu_AO
    m.ocean.τy .+= c.ρa * c.Cd_ao * ocn_frac .* vel_norm .* Δv_AO
    # Update ocean heatflux
    m.ocean.hflx_factor .= Δt * c.k/(c.ρi*c.L) .* (m.ocean.temp .- m.atmos.temp)
end