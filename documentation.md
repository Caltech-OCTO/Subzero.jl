# How to Create a Simulation: 

You have a lot of flexibility in designing a simulation within Subzero. First, you need to create all the pieces of your model, including the grid, the ocean, the atmosphere, the domain, and the floes. After than, you can use your model to create a simulation, where you will also specify physical constants and other runtime parameters.

## Building a Model
### Model Calculations Type: 

Simulations are designed to run with either Float64 or Float32 calculations. However, as of now, only runs in Float64 have been confirmed and are supported. When you are creating a model and a simulation, you need to create all elements in the same type that you are trying to run the model in. At the beginning of each example file, you will see the line: 

```julia 
const FT = Float64 
``` 

And then FT is used as an argument in the creation of the model and simulation. If all elements of the model and simulation are Float64, the simulation will run all calculations with Float64. Eventually, if all elements are created with Float32, the simulation will run all calculations with Float32. 

### Grid: 
The first thing that you need to create in your grid. For now, the ocean and atmosphere must be on the same grid. The only grid for now is a RegularRectilinearGrid, where every grid cell is a rectangle of the same size. You can create a regular rectilinear grid one of two ways. Both methods require specifying the minimum and maximum x and y points in meters. 

The first method has you provide the size of each grid cell with ∆x and ∆y arguments. That can be seen here: 

```julia 
RegRectilinearGrid( 
  FT, 
  xbounds = (-1e5, 1e5), 
  Ybounds = (0.0, 1e5) 
  ∆x = 2e4, 
  ∆y = 1e4, 
) 
``` 

Which will create a grid that has a width of 2e5m from –1e5m to 1e5m and each grid cell was a width of 2e4 and a height of 1e4. This would then create a 100 cell grid, 10 cells across and 10 cells wide. Note that if you provide grid cells dimensions that don’t divide the grid width and height evenly, the grid will be truncated to the nearest grid cell. 

The other way to create this same grid is to use: 

```julia 
RegRectilinearGrid(-1e5, 1e5, 0.0, 1e5, (10, 10)) 
``` 
by specifying the desired grid dimensions. Here you are guaranteed to always get a 10 by 10 grid, but do not have as fine-grained control over the size of your cells.  

A grid object has the following fields: dims, xg, yg, xc, and yc. The dims field is the number of rows and columns within the grid as defined by each grid cell. Therefore, dim = (number of cells in the y direction, number of cells in the x direction) as the y-direction cell count represents the number of rows in the grid, while the x-direction cell count is the number of columns. 

Additionally, RegRectilinearGrid is a concrete subtype of AbstractGrid. More concrete subtypes may be added in the future.  

### Ocean: 
The ocean here represents 2-dimensional vector fields. It includes the following: u velocity, v velocity, temperature, stress and forces in the x and y direction, the amount of area per cell covered in ice, and the hflx_factor per cell, which allows local heat flux calculations based on the difference in atmosphere and ocean temperature in each cell.  

### Atmosphere: 

### Domain: 
The domain includes both the boundaries of the simulation and the topographic features. It defines the areas where the floes cannot go. Before you can define a domain, you need to define each boundary wall. As of now, only rectangular boundaries around the edge of the grid are allowed so we need a north, east, south, and west boundary wall. We have four types of boundary wall: open, collision, periodic, and compression (not implemented yet).  

With an open wall, if a floe overlaps with the boundary wall at all, the floe is removed from the simulation. With a collision wall, floes collide with the wall, and it is calculated similarly to a collision with another floe, except that since the wall cannot break, or be moved, it has an idealized force factor. With a periodic wall, if a floe overlaps with the wall, a copy of that floe is created passing back into the domain through the opposite wall. These two floes are equivalent and are linked. We call the copy a “ghost floe.” So, for example, if a floe is pushed partially out of the domain through the south wall by a current, a copy of the floe will be created, re-entering through the north wall. If one wall in the domain is periodic, its opposite wall must also be periodic, that is north-south and east-west. Compression boundaries are walls that can move with a constant velocity towards the center of the domain. For example, if the west wall is a compression wall with a velocity of 0.1m/s it will move with a 0.1m/s u velocity and a 0m/s v velocity. These types of walls are for increasing pressure on the ice to investigate stress and strain on the floes.  

Here is an example of creating a set of boundary walls using the grid from above: 

```julia 
nboundary = PeriodicBoundary(grid, North()) 
sboundary = PeriodicBoundary(grid, South()) 
eboundary = CollisionBoundary(grid, East()) 
wboundary = OpenBoundary(grid, West()) 
``` 

Once we have defined our walls, we can create a domain as follows: 

```julia 
domain = Domain(nboundary, sboundary, eboundary, wboundary) 
``` 

However, if we want topography, we can also add it to the domain. To do that, we would create one, or more, TopographyElements. Let us consider a simple one, a square island: 

```julia 
island = [[[6e4, 4e4], [6e4, 4.5e4], [6.5e4, 4.5e4], [6.5e4, 4e4], [6e4, 4e4]]] 
topo_arr = StructVector([TopographyElement(t) for t in [island, topo1, topo2]]) 
``` 

We can then add an additional argument to the domain when it is created:  

```julia 
domain = Domain(nboundary, sboundary, eboundary, topo_arr) 
``` 

These are the two ways to define a domain in Subzero.  

For now, you will need to define ocean velocities that go around your topography, as Subzero does not change the given ocean velocities. This holds true on coupled scenarios as well. You need to make sure that if you do couple Subzero to Oceananigans that you set the ocean height to 0m at these locations so that the current floes around the topography. 

### Floes

### Making the Model
Once you have made all of the above components, you are now able to make a model. You will simply do that as follows:
```julia
model = Model(grid, ocean, atmos, domain, floe_arr)
```
This model object now holds all of the previously created components. You will get an error if you did not make all of your componenets with the same `FT` type. If you did make all objects with the same value of `FT`, either `Float64` or `Float32`, then the model now has that type.

## Building a Simulation
Your simulation will hold your model, as well as some runtime and physical parameters.

### Constants
Running the simulation requires known physical parameters, such as drag coefficents, densities, and the coriolis parameter. If you want to change any of the default values, you will need to create a `Constants` object. You can then keyword define any of the values that you want to change. All others will stay with their default value. For example, below we are creating a constants object that has a coefficent of friction set to 0, to run the simulation without friction. 

```julia
consts = Constants(μ = 0.0)
```
The list of constants and their default values are shown below:
| Constants     |    Meaning                         | Default Value |
| ------------- | ---------------------------------- | ------------- |
| ρi            | ice density                        | 920.0 kg/m^3  |
| ρo            | ocean densisty                     | 1027.0 kg/m^3 |
| ρa            | air densisty                       | 1.2 kg/m^3    |
| Cd_io         | ice-ocean drag coefficent          | 3e-3          |
| Cd_ia         | ice-air drag coefficent            | 1e-3          |
| Cd_ao         | air-ocean drag coefficent          | 1.25e-3       |
| f             | ocean coriolis frequency           | 1.4e-4 rad/s  |
| turnθ         | ocean turning angle                | 15π/180 rad   |
| L             | latent heat of freezing            | 2.93e5 J/kg   |
| k             | thermal conductivity of surface ice| 2.14 W/(mK)   |
| ν             | Poisson's ratio                    | 0.3           |
| μ             | Coefficent of friction             | 0.2           |
| E             | Young's Modulus                    | 6e6 N/m2      |

In particular, Young's Modulus is calculated usually calculated using the total floe area after floe initialization in the original Subzero code:
```julia
E = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
```

### Physical Process Settings 

### Timesteps
You have the ability to set the simulation's timestep in seconds using `∆t` and set the total number of timsteps the simulation will run for, `n∆t`. The default is `∆t = 10` seconds and `n∆t = 7500` timesteps. 

### Output Writers

### Reproducibility

### Creating the Simulation
