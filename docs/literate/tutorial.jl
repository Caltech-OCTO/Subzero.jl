# # Tutorial

# This is the user tutorial for Subzero.jl. It will walk you through the typical workflow of 
# buidling a discrete-element model (DEM) with Subzero.jl as well as running and plotting your
# simulatuion.


# ## Tutorial - copy-pasteable version

# ## Core steps of an Subzero.jl simulation


# The very first step of running a Subzero simulation is to bring the package into scope.
using Subzero  # bring package into scope

# ### Creating a grid

# Each Subzero model requires a grid object. The grid object defines the grid for the ocean
# and atmosphere. Ocean and atmosphere vectors (like `u` and `v` velocities) and tracers
# (like temperature) are recorded on these grid points and grid lines. The grid points are
# then used for interpolation for coupling between the ice, ocean, and atmosphere.

# All Subzero grid objects are concrete types of the abstract type (`AbstractGrid`)[@ref].
# Currently the only implemented concrete type is a (`RegRectilinearGrid`)[@ref]. If you are
# interested in implementing another type of grid, see the (developer documentation)["devdocs.md"].

# Here, we will go ahead and create an instance of `RegRectilinearGrid`. We need to specify
# the grid endpoints and either the number of grid cells in both directions, or the size of
# the grid cells. Here we will specity the number of grid cells in the x-direction, `Nx`, and
# in the y-direction, `Ny`.

grid = RegRectilinearGrid(; x0 = 0.0, xf = 1e5, y0 = -1e5, yf = 1e5, Nx = 10, Ny = 40)


# ### Creating a domain