"""
Plotting functions for Subzero Simulation
"""

function plot_domain!(plt, domain::CircleDomain)
    Θ = LinRange(0, 2π, 500)
    plot!(domain.centroid[1] .+ domain.radius*sin.(Θ),
          domain.centroid[2] .+ domain.radius*cos.(Θ),
          seriestype = [:shape], linecolor = :black, fillalpha = 0.0, lw = 2, legend=false)
end

function plot_domain!(plt, domain::RectangleDomain)
    north = domain.north.val
    south = domain.south.val
    east = domain.east.val
    west = domain.west.val
    plot!([west, west, east, east, west], [north, south, south, north, north],
           seriestype = [:shape], linecolor = :black, fillalpha = 0.0, lw = 2, legend=false)
end

function plot_sim(model, plt, time, PERIODIC) # Add Periodic
    # Setting up the Plots
    xgrid, ygrid = Subzero.grids_from_lines(model.grid.xc, model.grid.xc)
    quiver!(plt, vec(xgrid), vec(ygrid),
           quiver=(vec(model.ocean.u), vec(model.ocean.v)))
    plot_domain!(plt, model.domain)
    floe_coords = model.floes.coords
    floe_centroids = model.floes.centroid
    plot!([LG.Polygon(Subzero.translate(floe_coords[i], floe_centroids[i]))
           for i in eachindex(floe_coords)], fill = :lightblue)
    
    topos_coords = model.topos.coords
    topo_centroids = model.topos.centroid
    plot!([LG.Polygon(Subzero.translate(topos_coords[i], topo_centroids[i]))
           for i in eachindex(topos_coords)], fill = :grey)
    
    # I am not sure this is the most efficent way to do this
end