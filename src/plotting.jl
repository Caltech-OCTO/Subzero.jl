"""
Plotting functions for Subzero Simulation
"""

"""
grids_from_lines(xlines, ylines)

Creates x-grid and y-grid. Assume xlines has length n and ylines has length m. xgrid is the grid's xline vector repeated m times as rows in a mxn array and ygrid is the yline vector repeated n times as columns in a mxn vector. xlines and ylines are typically either xg and yg of xc and yc.
"""
function grids_from_lines(xlines, ylines)
    xgrid = repeat(reshape(xlines, 1, :), inner=(length(ylines),1))
    ygrid = repeat(ylines, outer = (1, length(xlines)))
    return xgrid, ygrid
end

"""
    domain_xycoords(domain::CircleDomain)

X-y coordinates for a rough circle in shape of given domain to be used for plotting.
Inputs:
        domain <CircleDomain>
Outputs:
        x-coordinates <Float Vector>
        y-coordinates <Float Vector>
"""
function domain_xycoords(domain::CircleDomain)
    Θ = LinRange(0, 2π, 500)
    return domain.centroid[1] .+ domain.radius*sin.(Θ),
           domain.centroid[2] .+ domain.radius*cos.(Θ)
end

"""
    domain_xycoords(domain::RectangleDomain)

X-y coordinates for a rectangle in shape of given domain to be used for plotting.
Inputs:
        domain <RectangleDomain>
Outputs:
        x-coordinates <Float Vector>
        y-coordinates <Float Vector>
"""
function domain_xycoords(domain::RectangleDomain)
    north = domain.north.val
    south = domain.south.val
    east = domain.east.val
    west = domain.west.val
    return [west, west, east, east, west], [north, south, south, north, north]
end

"""
    setup_plot(model::Model)

Plots grid lines, model domain, and model topology on a plot with a ratio determined by the grid coordinates. This plot will be used throughout simulation so these qualities only need to be plotted/set once as they cannot change during the simulation.

Inputs:
        model <Model>
Output:
        plt <Plots.Plot> with gridlines and aspect ratio set and model domain and topology plotted
"""
function setup_plot(model::Model)
    # Plot Grid Lines and set image ratio
    xmax = model.grid.xg[end]
    ymax = model.grid.yg[end]
    ratio = ymax/xmax
    plt = Plots.plot(xlims = (model.grid.xg[1], xmax),
                     ylims = (model.grid.yg[1], ymax), size = (1200, 1200), aspect_ratio=ratio)
    # Plot Domain Border
    domainx, domainy = domain_xycoords(model.domain)
    plot!(plt, domainx, domainy, seriestype = [:shape], linecolor = :black, 
            fillalpha = 0.0, lw = 2, legend=false)

    # Plot Topology
    topos_coords = model.topos.coords
    topo_centroids = model.topos.centroid
    plot!(plt, [LG.Polygon(Subzero.translate(topos_coords[i],
            topo_centroids[i])) for i in eachindex(topos_coords)], fill = :grey)
    return plt
end

"""
    plot_sim(model, plt)

Plot ocean velocities and floes for current model time step using given plot.
    
Inputs:
        model   <Model>
        plt     <Plots.plt>
        time    <Int> model timestep
Outpits:
        Save new plot
"""
function plot_sim(model, plt, time)
    # Plot Ocean Vector Field - also clears previous plot
    xgrid, ygrid = Subzero.grids_from_lines(model.grid.xc, model.grid.xc)
    plt_new = quiver(plt, vec(xgrid), vec(ygrid),
            quiver=(vec(model.ocean.u), vec(model.ocean.v)), color = :lightgrey)

    # Plot Floes --> only plot "alive" floes
    floe_coords = model.floes.coords
    floe_centroids = model.floes.centroid
    floe_alive = model.floes.alive
    plot!(plt_new, [LG.Polygon(Subzero.translate(floe_coords[i],
          floe_centroids[i])) for i in eachindex(floe_coords) if floe_alive[i]],
          fill = :lightblue)
          
    # Save plot
    Plots.savefig(plt_new, "figs/plot_$time.png")
end