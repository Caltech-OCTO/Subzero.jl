"""
Functions needed for coupling the ice, ocean, and atmosphere.
"""

#-------------- Monte Carlo Point Calculations --------------#

"""
    find_cell_indices(xp, yp, grid::RegRectilinearGrid)

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
function find_cell_indices(xp, yp, grid::RegRectilinearGrid)
    xidx = floor.(Int, (xp .- grid.xg[1])/(grid.xg[2] - grid.xg[1]) .+ 0.5) .+ 1
    yidx = floor.(Int, (yp .- grid.yg[1])/(grid.yg[2] - grid.yg[1]) .+ 0.5) .+ 1
    return xidx, yidx
end

function find_cell_index(xp, yp, grid::RegRectilinearGrid)
    xidx = floor(Int, (xp - grid.xg[1])/(grid.xg[2] - grid.xg[1]) + 0.5) + 1
    yidx = floor(Int, (yp - grid.yg[1])/(grid.yg[2] - grid.yg[1]) + 0.5) + 1
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

function in_bounds(
    xr,
    yr,
    grid,
    ::NonPeriodicBoundary,
    ::NonPeriodicBoundary,
)
    return (grid.xg[1] <= xr <= grid.xg[end]) && 
        (grid.yg[1] <= yr <= grid.yg[end])
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

function in_bounds(
    xr,
    yr,
    grid,
    ::NonPeriodicBoundary,
    ::PeriodicBoundary,
)
    return grid.xg[1] <= xr <= grid.xg[end]
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

function in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::NonPeriodicBoundary,
)
    return grid.yg[1] <= yr <= grid.yg[end]
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

function in_bounds(
    xr,
    yr,
    grid,
    ::PeriodicBoundary,
    ::PeriodicBoundary,
)
    return true
end

"""
    calc_mc_values(
        floe,
        grid,
        domain,
    )

Calculates monte carlo point's cartesian coordiantes, polar coordiantes,
velocity and index within the grid. 
Inputs:
    floe    <Floe> single floe
    grid    <AbstractGrid> model grid
    domain  <Domain> model domain
Outputs:
    mc_cart     <Matrix{AbstractFloat}> cartesian coordinates for model
                    coordinates - nx2 matrix of monte carlo coordinates where
                    first column is the x-coords and the second column is the
                    y-coords
    mc_polar    <Matrix{AbstractFloat}> polar coordinates for model coordinates
                    - nx2 matrix of monte carlo coordinates where first column
                    is the radius and the second column is the angle
    mc_vel     <Matrix{AbstractFloat}> velocity at each monte carlo point - nx2
                    matrix where the first column is the u-velocity and the
                    second is the v-velocity for n monte carlo points
    mc_grid_idx <Matrix{Int}> index of monte carlo points within the grid - nx2
                    matrix of indices where the first column is the grid column
                    index and the second column is the grid row index
"""
function calc_mc_values!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    grid,
    domain,
    mc_cart,
    mc_polar,
    mc_vel,
    mc_grid_idx,
) where {FT}
    # Translate/rotate monte carlo points to floe location/orientation
    α = floe.α
    j = 0  # index in output array
    for i in eachindex(floe.mc_x)  # index of monte carlo points
        mc_px = cos(α)*floe.mc_x[i] - sin(α)*floe.mc_y[i]  # at origin
        mc_py = sin(α)*floe.mc_x[i] + cos(α)*floe.mc_y[i]  # at origin
        mc_x = mc_px + floe.centroid[1]  # at centroid
        mc_y = mc_py + floe.centroid[2]  # at centroid
        # If point is in bounds, continue to find rest of values
        if in_bounds(mc_x, mc_y, grid, domain.east, domain.north)
            j += 1  # if added to outputs, move to next index in output array
            mc_cart[j, 1] = mc_x
            mc_cart[j, 2] = mc_y
            mc_polar[j, 1] = sqrt(mc_px^2 + mc_py^2)
            mc_polar[j, 2] = atan(mc_py, mc_px)
            mc_vel[j, 1] = floe.u - floe.ξ * mc_polar[j, 1] * sin(mc_polar[j, 2])
            mc_vel[j, 2] = floe.v + floe.ξ * mc_polar[j, 1] * cos(mc_polar[j, 2])
            mc_grid_idx[j, 1], mc_grid_idx[j, 2] = find_cell_index(
                mc_cart[j, 1],
                mc_cart[j, 2],
                grid,
            )
        end
    end
    return j  # last spot in output arrays that has values for this floe
end

#-------------- Interpolation of Ocean and Atmosphere --------------#

"""
    find_interp_knots(
        point_idx,
        glines,
        Δd::Int,
        ::PeriodicBoundary,
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
    glines::Vector{FT},
    Δd::Int,
    ::PeriodicBoundary,
) where {FT<:AbstractFloat} 
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
    glines::Vector{FT},
    Δd::Int,
    ::NonPeriodicBoundary,
) where {FT}
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
function mc_interpolation!(
    npoints,
    mc_cart,
    mc_grid_idx,
    grid,
    domain,
    atmos,
    ocean,
    coupling_settings,
    uatm,
    vatm,
    uocn,
    vocn,
    hflx_factor,
)
    mc_xr = @view mc_cart[1:npoints, 1]
    mc_yr = @view mc_cart[1:npoints, 2]
    mc_cols = @view mc_grid_idx[1:npoints, 1]
    mc_rows = @view mc_grid_idx[1:npoints, 2]

    # Find knots and indices of knots for monte carlo interpolation
    xknots, xknot_idx = find_interp_knots(
        mc_cols,
        grid.xg,
        coupling_settings.Δd,
        domain.east,
    )
    yknots, yknot_idx = find_interp_knots(
        mc_rows,
        grid.yg,
        coupling_settings.Δd,
        domain.north,
    )
    knots = (yknots, xknots)
    # Atmos Interpolation for Monte Carlo Points
    uatm_interp = linear_interpolation(
        knots,
        @view(atmos.u[yknot_idx, xknot_idx]),
    )
    
    vatm_interp = linear_interpolation(
        knots,
        @view(atmos.v[yknot_idx, xknot_idx]),
    )
    
    # Ocean Interpolation for Monte Carlo Points
    uocn_interp = linear_interpolation(
        knots,
        @view(ocean.u[yknot_idx, xknot_idx]),
    )

    vocn_interp = linear_interpolation(
        knots,
        @view(ocean.v[yknot_idx, xknot_idx]),
    )

    # Heatflux Factor Interpolatoin for Monte Carlo Points
    hflx_interp = linear_interpolation(
        knots,
        @view(ocean.hflx_factor[yknot_idx, xknot_idx]),
    )

    for i in eachindex(mc_xr)
        uatm[i] = uatm_interp(mc_yr[i], mc_xr[i]) 
        vatm[i] = vatm_interp(mc_yr[i], mc_xr[i])
        uocn[i] = uocn_interp(mc_yr[i], mc_xr[i])
        vocn[i] = vocn_interp(mc_yr[i], mc_xr[i])
        hflx_factor[i] = hflx_interp(mc_yr[i], mc_xr[i])
    end
    return
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

function add_point!(scell::OceanStressCell, floeidx, τx, τy, Δx, Δy)
    if isempty(scell.floeidx) || scell.floeidx[end] != floeidx
        push!(scell.floeidx, floeidx)
        push!(scell.τx, τx)
        push!(scell.τy, τy)
        push!(scell.npoints, 1)
        push!(scell.trans_vec, SVector{2}([Δx, Δy]))
    else
        scell.τx[end] += τx
        scell.τy[end] += τy
        scell.npoints[end] += 1
    end

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
    )

Add force from the ice on ocean to ocean force fields (fx & fy) for each grid
cell and update ocean sea ice area fraction (si_area), representing total area
of sea ice in a given cell.
Inputs:
    mc_grid_idx <Matrix{Int}> index of monte carlo points within the grid - nx2
                    matrix of indices where the first column is the grid column
                    index and the second column is the grid row index
    τx_ocn      <Vector{Float}> x-stress caused by ocean on each corresponding
                    monte carlo point - negative is stress on the ocean
    τy_ocn      <Vector{Float}> y-stress caused by ocean on each corresponding
                    monte carlo point - negative is stress on the ocean
    floe        <Floe> single floe within model
    ocean       <Ocean> model's ocean
    grid        <AbstractGrid> model's grid
    ns_bound    <AbstractBoundary> north or south boundary of domain -
                    dispatches on periodic vs non-periodic
    ew_bound    <AbstractBoundary> east or west boundary of domain -
                    dispatches on periodic vs non-periodic
    spinlock    <Thread.SpinLock>
Outputs:
    None. Ocean fields updated in-place. Note that the needed additions take
    place within a lock. Since multiple floes can be within one grid cell, this
    critical section needs to be locked for when the code is run with multiple
    threads to prevent race conditions.
"""
function aggregate_grid_force!(
    floeidx,
    npoints,
    mc_grid_idx,
    τx_ocn::Vector{FT},
    τy_ocn,
    grid,
    ns_bound,
    ew_bound,
    scells,
) where {FT}
    mc_cols = @view mc_grid_idx[1:npoints, 1]
    mc_rows = @view mc_grid_idx[1:npoints, 2]
    # Use dictionary to sort stress values by grid cell
    for i in eachindex(mc_cols)
        row = shift_cell_idx(mc_rows[i], grid.dims[1] + 1, ns_bound)
        col = shift_cell_idx(mc_cols[i], grid.dims[2] + 1, ew_bound)
        Δx = (col - mc_cols[i]) * (grid.xg[2] - grid.xg[1])
        Δy = (row - mc_rows[i]) * (grid.yg[2] - grid.yg[1])
        add_point!(scells[row, col], floeidx, -τx_ocn[i], -τy_ocn[i], Δx, Δy)
    end
    return
end

function sum_grid_force!(
    floes::StructArray{Floe{FT}},
    grid::RegRectilinearGrid,
    atmos,
    ocean,
    consts,
    ns_bound,
    ew_bound,
) where {FT}
    # Determine force from floe on each grid cell it is in
    cell_area = (grid.xg[2] - grid.xg[1]) * (grid.yg[2] - grid.yg[1])
    Threads.@threads for cartidx in CartesianIndices(ocean.scells)
        ocean.τx[cartidx] = FT(0)
        ocean.τy[cartidx] = FT(0)
        ocean.si_frac[cartidx] = FT(0)
        τocn = ocean.scells[cartidx]
        if !isempty(τocn.floeidx)
            # Coordinates of grid cell
            cell_coords = center_cell_coords(
                cartidx[2],
                cartidx[1],
                grid,
                ns_bound,
                ew_bound
            )
            cell_poly = LG.Polygon(cell_coords)
            for i in eachindex(τocn.floeidx)
                floe_coords = translate(
                    floes.coords[τocn.floeidx[i]],
                    τocn.trans_vec[i][1],
                    τocn.trans_vec[i][2],
                )
                floe_poly = LG.Polygon(floe_coords)
                floe_area_in_cell = LG.area(LG.intersection(cell_poly, floe_poly))::FT
                if floe_area_in_cell > 0
                    # Add forces and area to ocean fields
                    ocean.τx[cartidx] += (τocn.τx[i]/τocn.npoints[i]) * floe_area_in_cell
                    ocean.τy[cartidx] += (τocn.τy[i]/τocn.npoints[i]) * floe_area_in_cell
                    ocean.si_frac[cartidx] += floe_area_in_cell
                end
            end
            if ocean.si_frac[cartidx] > 0
                # Divide by total floe area in cell to get ocean stress
                ocean.τx[cartidx] /= ocean.si_frac[cartidx]
                ocean.τy[cartidx] /= ocean.si_frac[cartidx]
                # Divide by cell area to get sea-ice fraction
                ocean.si_frac[cartidx] /= cell_area
            end
        end
        Δu_AO = atmos.u[cartidx] - ocean.u[cartidx]
        Δv_AO = atmos.v[cartidx] - ocean.v[cartidx]
        ocn_frac = 1 - ocean.si_frac[cartidx]
        norm = sqrt(Δu_AO^2 + Δv_AO^2)
        ocean.τx[cartidx] += consts.ρa * consts.Cd_ao * ocn_frac * norm * Δu_AO
        ocean.τy[cartidx] += consts.ρa * consts.Cd_ao * ocn_frac * norm * Δv_AO
        # Not sure this is where the heatflux should be??
        ocean.hflx_factor[cartidx] = Δt * consts.k/(consts.ρi*consts.L) *
            (ocean.temp[cartidx] - m.atmos.temp[cartidx])
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
function calc_ocean_forcing!(
    floe::Union{Floe{FT}, LazyRow{Floe{FT}}},
    mc_vel,
    uocn,
    vocn,
    c,  # constants
    τx_ocn,
    τy_ocn,
    τx_pressure∇,
    τy_pressure∇,
) where {FT}
    ma_ratio = floe.mass/floe.area
    # Stress on ice from ocean
    for i in eachindex(uocn)
        Δu_OI = uocn[i] - mc_vel[i, 1]
        Δv_OI = vocn[i] - mc_vel[i, 2]
        norm = sqrt(Δu_OI^2 + Δv_OI^2)
        τx_ocn[i] = c.ρo*c.Cd_io * norm * (cos(c.turnθ) * Δu_OI - sin(c.turnθ) * Δv_OI)
        τy_ocn[i] = c.ρo*c.Cd_io * norm * (sin(c.turnθ) * Δu_OI + cos(c.turnθ) * Δv_OI)
        τx_pressure∇[i] = -ma_ratio * c.f * vocn[i]
        τy_pressure∇[i] = ma_ratio * c.f * uocn[i]
    end
    return
end

"""
    timestep_coupling!(
        model,
        consts,
        coupling_settings,
        spinlock,
    )

Calculates the effects of the ocean and atmosphere on the ice and the effects of
the ice on the ocean.
Inputs:
    model               <Model> model
    consts              <Constants> constants used in simulation
    coupling_settings   <CouplingSettings> settings for coupling
    spinlock            <Thread.SpinLock>
Outputs:
    None. Updates each floe's ocean/atmosphere forcings (fxOA, fyOA, torqueOA)
    and calculates stresses on each ocean grid cell.          
"""
function timestep_coupling!(
    floes::StructArray{Floe{FT}},
    grid,
    domain,
    ocean,
    atmos,
    consts,
    coupling_settings,
) where {FT<:AbstractFloat}
    # Allocations
    # Monte carlo point values
    max_n_mc = maximum(length, floes.mc_x)
    mc_cart = Matrix{FT}(undef, max_n_mc, 2)
    mc_polar = Matrix{FT}(undef, max_n_mc, 2)
    mc_vel = Matrix{FT}(undef, max_n_mc, 2)
    mc_grid_idx = Matrix{Int}(undef, max_n_mc, 2)
    # Ocean/atmosphere values
    uatm = Vector{FT}(undef, max_n_mc)
    vatm = Vector{FT}(undef, max_n_mc)
    uocn = Vector{FT}(undef, max_n_mc)
    vocn = Vector{FT}(undef, max_n_mc)
    hflx_factor = Vector{FT}(undef, max_n_mc)
    # Forces and Torques
    τx_ocn = Vector{FT}(undef, max_n_mc)
    τy_ocn = Vector{FT}(undef, max_n_mc)
    τx_pressure∇ = Vector{FT}(undef, max_n_mc)
    τy_pressure∇ = Vector{FT}(undef, max_n_mc)
    τx = Vector{FT}(undef, max_n_mc)
    τy = Vector{FT}(undef, max_n_mc)
    τtrq = Vector{FT}(undef, max_n_mc)
    # Clear ocean forces and area fractions
    empty!.(ocean.scells)

    # Calcualte coupling for each floe
    for i in eachindex(floes)
        # Find monte carlo point peroperties
        npoints = calc_mc_values!(
            LazyRow(floes, i),
            grid,
            domain,
            mc_cart,
            mc_polar,
            mc_vel,
            mc_grid_idx,
        )

        # Interpolate ocean and atmosphere values onto monte carlo points
        mc_interpolation!(
            npoints,
            mc_cart,
            mc_grid_idx,
            grid,
            domain,
            atmos,
            ocean,
            coupling_settings,
            uatm,
            vatm,
            uocn,
            vocn,
            hflx_factor,
        )
        # Calculate effects on the floe
        floes.hflx_factor[i] = mean(@view hflx_factor[1:npoints])
        τx_atm, τy_atm = @views calc_atmosphere_forcing(
            uatm[1:npoints],
            vatm[1:npoints],
            consts,
        )
        @views calc_ocean_forcing!(
            LazyRow(floes, i),
            mc_vel[1:npoints, :],
            uocn[1:npoints],
            vocn[1:npoints],
            consts,
            τx_ocn,
            τy_ocn,
            τx_pressure∇,
            τy_pressure∇,
        )

        # Update ocean with force from floe per grid cell
        if coupling_settings.calc_ocnτ_on
            # Calculate effects of floe on ocean
            aggregate_grid_force!(
                i,
                npoints,
                mc_grid_idx,
                τx_ocn,
                τy_ocn,
                grid,
                domain.north,
                domain.east,
                ocean.scells,
            )
        end

        # Sum stresses and find torques
        xcoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.v[i]
        ycoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.u[i]
        for j in 1:npoints
            τx[j] = τx_atm + τx_pressure∇[j] + τx_ocn[j]
            τy[j] = τy_atm + τy_pressure∇[j] + τy_ocn[j]
            # Calculate torque
            τtrq[j] = (-τx[j] * sin(mc_polar[j, 2]) + τy[j] * cos(mc_polar[j, 2])) * mc_polar[j, 1]
            # Add coriolis force to total foces - does not contribute to torque
            τx[j] += xcoriolis
            τy[j] -= ycoriolis
        end

        # Average forces on ice floe
        floes.fxOA[i] = mean(@view τx[1:npoints]) * floes.area[i]
        floes.fyOA[i] = mean(@view τy[1:npoints]) * floes.area[i]
        floes.trqOA[i] = mean(@view τtrq[1:npoints]) * floes.area[i]
    end
    if coupling_settings.calc_ocnτ_on
        sum_grid_force!(
            floes,
            grid,
            atmos,
            ocean,
            consts,
            domain.north,
            domain.east,
        )
    end
end
