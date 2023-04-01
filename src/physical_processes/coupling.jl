"""
Functions needed for coupling the ice, ocean, and atmosphere.
"""

#-------------- Monte Carlo Point Calculations --------------#

"""
    find_cell_index(xp, yp, grid::RegRectilinearGrid)

Find index of the cell centered on grid lines of the given RegRectilinearGrid
that the given x-coordinate and y-coordinate fall within.
This cell is centered around the grid lines, so it is a shifted grid cell
by half a cell. Method depends on grid being a regular rectilinear grid.
Inputs:
    xp      <AbstractFloat> x-coordinates of point
    yp      <AbstractFloat> y-coordinate of point
    grid    <RegRectilinearGrid> simulation grid
Outputs:
    xidx    <AbstractFloat> x-index of grid cell (cented on grid lines) x-point
                is within
    yidx    <AbstractFloat> y-index of grid cell (cented on grid lines) y-point
                is within
Note:
    Points can be outside of the grid, so index can be less than 1 or greater
    than the number of grid lines in a given direction.
"""
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
    mc_grid_idx,
    grid,
    domain,
    atmos,
    ocean,
    coupling_settings,
)
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

    return uatm_interp, vatm_interp, uocn_interp, vocn_interp, hflx_interp
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
    mc_xr, 
    mc_yr,
    upoint,
    vpoint,
    uatm_interp,
    vatm_interp,
    c,  # constants
)
    # Atmosphere velocities at monte carlo point
    uatm = uatm_interp(mc_yr, mc_xr) 
    vatm = vatm_interp(mc_yr, mc_xr)

    # Stress on ice from atmopshere
    Δu_AI = uatm - upoint
    Δv_AI = vatm - vpoint
    norm = sqrt(Δu_AI^2 + Δv_AI^2)
    τx_atm = c.ρa * c.Cd_ia * norm * Δu_AI
    τy_atm = c.ρa * c.Cd_ia * norm * Δv_AI
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
    mc_xr,
    mc_yr,
    upoint,
    vpoint,
    uocn_interp,
    vocn_interp,
    hflx_interp,
    ma_ratio,
    c,  # constants
)
    uocn = uocn_interp(mc_yr, mc_xr)
    vocn = vocn_interp(mc_yr, mc_xr)
    hflx_factor = hflx_interp(mc_yr, mc_xr)
    Δu_OI = uocn - upoint
    Δv_OI = vocn - vpoint
    norm = sqrt(Δu_OI^2 + Δv_OI^2)
    τx_ocn = c.ρo*c.Cd_io * norm * (cos(c.turnθ) * Δu_OI - sin(c.turnθ) * Δv_OI)
    τy_ocn = c.ρo*c.Cd_io * norm * (sin(c.turnθ) * Δu_OI + cos(c.turnθ) * Δv_OI)
    τx_pressure∇ = -ma_ratio * c.f * vocn
    τy_pressure∇ = ma_ratio * c.f * uocn
    return τx_ocn, τy_ocn, τx_pressure∇, τy_pressure∇, hflx_factor
end

function calc_one_way_coupling!(
    floes::StructArray{Floe{FT}},
    grid,
    atmos,
    ocean,
    domain,
    coupling_settings,
    consts,
) where {FT}
    max_n_mc = maximum(length, floes.mc_x)
    mc_cart = Matrix{FT}(undef, max_n_mc, 2)
    mc_grid_idx = Matrix{Int}(undef, max_n_mc, 2)
    for i in eachindex(floes)
        # Monte carlo point cartesian coordinates and grid cell indices
        npoints = calc_mc_values!(
            LazyRow(floes, i),
            grid,
            domain,
            mc_cart,
            mc_grid_idx,
        )
        # Interpolaters for ocean and atmosphere
        uatm_int, vatm_int, uocn_int, vocn_int, hflx_int = mc_interpolation!(
            npoints,
            mc_grid_idx,
            grid,
            domain,
            atmos,
            ocean,
            coupling_settings,
        )

        # Add coriolis stress to total stress - same for every monte carlo point
        xcoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.v[i]
        ycoriolis = (floes.mass[i]/floes.area[i]) * consts.f * floes.u[i]
        tot_τx = npoints * xcoriolis
        tot_τy = -npoints * ycoriolis
        tot_τtrq = FT(0)
        tot_hflx_factor = FT(0)
        ma_ratio = floes.mass[i]/floes.area[i]
        for j in 1:npoints
            # Monte carlo point properties
            xcentered = mc_cart[j, 1] - floes.centroid[i][1]

            ycentered = mc_cart[j, 2] - floes.centroid[i][2]
            θ = atan(ycentered, xcentered)
            rad = sqrt(xcentered^2 + ycentered^2)
            upoint = floes.u[i] - floes.ξ[i] * rad * sin(θ)
            vpoint = floes.v[i] + floes.ξ[i] * rad * cos(θ)
            # Stress at monte carlo point from ocean and atmosphere
            τx_atm, τy_atm = calc_atmosphere_forcing(
                mc_cart[j, 1], 
                mc_cart[j, 2],
                upoint,
                vpoint,
                uatm_int,
                vatm_int,
                consts,
            )
            τx_ocn, τy_ocn, τx_p∇, τy_p∇, hflx_factor = calc_ocean_forcing!(
                mc_cart[j, 1],
                mc_cart[j, 2],
                upoint,
                vpoint,
                uocn_int,
                vocn_int,
                hflx_int,
                ma_ratio,
                consts,
            )
            τx = τx_atm + τx_p∇ + τx_ocn
            τy = τy_atm + τy_p∇ + τy_ocn
            # Torque at monte carlo point
            τtrq = (-τx * sin(θ) + τy * cos(θ)) * rad
            # Add values to total stresses
            tot_τx += τx
            tot_τy += τy
            tot_τtrq += τtrq
            tot_hflx_factor += hflx_factor

            floe_to_grid_info!(
                i,
                mc_grid_idx[j, 2],  # row
                mc_grid_idx[j, 1],  # column
                τx_ocn,
                τy_ocn,
                grid,
                domain,
                ocean.scells,
            )
        end
        # Average forces on ice floe
        floes.fxOA[i] = tot_τx/npoints * floes.area[i]
        floes.fyOA[i] = tot_τy/npoints * floes.area[i]
        floes.trqOA[i] = tot_τtrq/npoints * floes.area[i]
        floes.hflx_factor[i] = tot_hflx_factor/npoints
    end
end



"""
    floe_to_grid_info!(
        floeidx,
        row,
        col,
        τx_ocn::FT,
        τy_ocn::FT,
        grid,
        domain,
        scells,
    )

Add force from the ice on ocean to ocean force fields (fx & fy) for each grid
cell and update ocean sea ice area fraction (si_area), representing total area
of sea ice in a given cell.
Inputs:
    floeidx     <Int> index of floe within model's floe array
    row         <Int> 
    col         <Int>
    τx_ocn      <AbstractFloat> x-stress caused by ocean on monte carlo point
    τy_ocn      <AbstractFloat> y-stress caused by ocean on monte carlo point

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
function floe_to_grid_info!(
    floeidx,
    row,
    col,
    τx_ocn::FT,
    τy_ocn::FT,
    grid,
    domain,
    scells,
) where {FT}
    # Use dictionary to sort stress values by grid cell
    shifted_row = shift_cell_idx(row, grid.dims[1] + 1, domain.north)
    shifted_col = shift_cell_idx(col, grid.dims[2] + 1, domain.east)
    Δx = (shifted_col - col) * (grid.xg[2] - grid.xg[1])
    Δy = (shifted_row - row) * (grid.yg[2] - grid.yg[1])
    add_point!(
        scells[shifted_row, shifted_col],
        floeidx,
        -τx_ocn,
        -τy_ocn,
        Δx,
        Δy,
    )
    return
end

"""
    calc_two_way_coupling!(
        floes::StructArray{Floe{FT}},
        grid::RegRectilinearGrid,
        atmos,
        ocean,
        domain,
        consts,
        Δt,
    )

Calculate effects of ice and atmosphere on the ocean and update ocean stress
fields and sea ice fraction.
Inputs:
    floes       <StructArray{Floe}> model's floes
    ocean       <Ocean> model's ocean
    grid        <AbstractGrid> model's grid
    domain      <Domain> model's domain
    consts      <Constants> model's constantd
    Δt          <Int> simulation's timestep in seconds
Output:
    None. Update's ocean's stress fields and heatflux factor field. 
"""
function calc_two_way_coupling!(
    floes::StructArray{Floe{FT}},
    grid::RegRectilinearGrid,
    atmos,
    ocean,
    domain,
    consts,
    Δt,
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
                domain.north,
                domain.east
            )
            cell_poly = LG.Polygon(cell_coords)
            for i in eachindex(τocn.floeidx)
                floe_coords = translate(
                    floes.coords[τocn.floeidx[i]],
                    τocn.trans_vec[i][1],
                    τocn.trans_vec[i][2],
                )
                floe_poly = LG.Polygon(floe_coords)
                floe_area_in_cell = LG.area(LG.intersection(
                    cell_poly,
                    floe_poly,
                ))::FT
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
            (ocean.temp[cartidx] - atmos.temp[cartidx])
    end
    return
end


"""
    timestep_coupling!(
        model,
        Δt,
        consts,
        coupling_settings,
    )

Calculates the effects of the ocean and atmosphere on the ice and the effects of
the ice and atmosphere on the ocean if the coupling is two-way.
Inputs:
    model               <Model> model
    Δt                  <Int> length of timestep in seconds
    consts              <Constants> constants used in simulation
    coupling_settings   <CouplingSettings> settings for coupling
Outputs:
    None. Updates each floe's ocean/atmosphere forcings (fxOA, fyOA, torqueOA)
    and calculates stresses on each ocean grid cell from ice and atmosphere if
    two-way coupling is enabled in coupling_settings       
"""
function timestep_coupling!(
    model,
    Δt,
    consts,
    coupling_settings,
)
    calc_one_way_coupling!(
        model.floes,
        model.grid,
        model.atmos,
        model.ocean,
        model.domain,
        coupling_settings,
        consts,
    )
    if coupling_settings.calc_ocnτ_on
        calc_two_way_coupling!(
            model.floes,
            model.grid,
            model.atmos,
            model.ocean,
            model.domain,
            consts,
            Δt,
        )
    end
    return
end
