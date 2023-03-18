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

Floes are quite complex objects as they have a lot of needed fields. Here we will talk about a floe struct's fields, as well as how to create a configuration of floes to start your simulation.

One important thing to know is that a floe's coordinates are represented by a `PolyVec`, which is a shorthand for a vector of a vector of a vector of floats. This sounds complicated, but it is simply a way of representing a polygon's coordinates. A Polygon's coordinates are of the form below, where the xy-coordinates are the exterior border of the floe and the wz-coordinates, or any other following sets of coordinates, describe holes within the floe:

```julia
coords = [
  [[x1, y1], [x2, y2], ..., [xn, yn], [x1, y1]],  # Exterior vertices of the polygon represented as a list of cartesian points
  [[w1, z1], [w2, z2], ..., [wn, zn], [w1, z1]],  # Interior holes of the polygon represented as a list of cartesian points
  ...,  # Additional holes within the polygon represented as a list of cartesian points
 ]
 ```
 We will use the term `PolyVec` to describe this form of coordiantes and you will see it in the code if you take a look at the code base. It is also the form that floe coordinates are saved in output files.

It is recomeneded that you use the `initialize_floe_field` to create your starting configuration on floes. There are two ways to use this function. One is to provide a list of `PolyVecs` representing a list of the coordinates of all of the floes you want in the initial state. This will initialize all of the given polygons specified as floes. The other is to provide a number of floes and a concentration over a specific area. This will create a starting floe field using Voronoi Tesselation that aims to achieve the requested number of floes and concentrations.

Both of these functions share most arguments. They are as follows (with default values if they exist):
- domain, which is the model's domain so that the floes can be fit into the open space using polygon intersections/differences
- hmean, which is the mean height of all floes created
- Δh, which is the maximum potential height difference from `hmean` between floes
- min_floe_area = 0.0, which is the minimum floe area for any floes initialized and smaller floes will not be added to the list of initial floes
- ρi = 920.0, which is the density of ice in kg/m^3
- mc_n::Int = 1000, which is the number of monte carlo points desired for each floe
- nhistory::Int = 1000, which is the length of stress history for each floe
- rng = Xoshiro(), which is a random number generator so that floe creation is reproducible if a seeded random number generator is provided

Note that all arguments with default values are optional keyword arguments.

Here is an example of creating a small floe field using the version of `initialize_floe_field` that takes in lists of `PolyVec`s. For a real simulation, you would probably generate a list of coordinates and read them in from a file but we create these by hand here for simplicity. Assume we have already created a `domain`.
```julia
floe1 = [[[6e4, 2e4], [6e4, 5e4], [9e4, 5e4], [9e4, 2e4], [6e4, 2e4]]]
floe2 = [[[5.5e4, 2e4], [5.25e4, 4e4], [5.75e4, 4e4], [5.5e4, 2e4]]]
floe_field = initialize_floe_field(
  [floe1, floe2],
  domain,
  0.25,  # mean height of 0.25
  0.0,  # all floes will be the same height
  rng = Xoshiro(1),
  nhistory = 100,
)
```

Now here is an an example of creating a large floe field using the version of `initialize_floe_field` that uses Voronoi tesselation. Again assume we have already created a `domain`.
```julia
floe_arr = initialize_floe_field(
    100,  # attempt to initialize 100 floes
    [1.0; 0.0],  # the top half of the domain is fully packed and the bottom has no floes
    domain,
    0.25,  # mean height of 0.25
    0.10,  # floe heights will range from 0.15-0.35
    min_floe_area = 1e7,
    rng = Xoshiro(1),
    nhistory = 1000,
)
```
We now focus on the first two arguments. The first is the number of floes to attempt to create with Voronoi tesselation. We are not guarenteed to get exactly that number. It depends on the amount of open space in the domain and the generation of random seed points. For example, if the domain is filled with lots of topography and islands, it will be more difficult to hit the exact number of floes requested. However, it will be in the ballpark. The second argument is the concentrations, which is a matrix. We can split the domain into quadrents that are the same shape at matrix and then request concentrations of ice in each of those quadrents equal to the corresponding value in the concentrations matrix. The other arguments are the same as in the floe coordinate version on the function.

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
Subzero allows you to turn on and off various physical processes, as well as change various settings for these physical processes. The physical processes availible are: coupling, collisions, fractures, and simplfication. We will be adding corner fracturing, packing, rafting, ridging, and welding. For each of these physical processes, you can create a settings object, and change the parameters within that object to change the physical process to meet your simulation needs. 
#### Coupling Settings
`CouplingSettings` changes how Subzero uses the ocean and atmosphere two-dimensional vector fields within the `Ocean` and `Atmos` structs, whether those values are static user-provided values of values that are updated every timesteps through coulpling with Oceananigans. The availible settings are:
- coupling_on, which turns this process on and off
- ∆t, which sets the number of timesteps between this process running
- ∆d, which sets the number of ocean/atmosphere grid cells around a floe to consider when interpolating ocean/atmosphere forcings
- mc_n, which sets the number of monte carlo points for each floe that are used for the interpolation
- calc_ocnτ_on, which turns on and off the calculation of the ice effects on the ocean and stores output in the ocean's fx, fy, and si_area fields

Here is an example of creating your own coupling settings, using the default values:
```julia
couple_settings = CouplingSettings(
  coupling_on = true,
  Δt = 10,
  Δd = 1,
  mc_n = 1000,
  calc_ocnτ_on = false,
)
Note that you only need to provide values to the fields that you wish to change from the defaults.
```
#### Collision Settings
`CollisionSettings` changes Subzero's floe-floe and floe-domain collisions. The availible settings are:
- collisions_on, which turns this process on and off
- floe_floe_max_overlap, which sets the maximum fraction of floe-floe overlap before fusing the two floes together 
- floe_domain_max_overlap, which sets the maximum fraction for floe-domain overlap before removing the floe from the simulation

Here is an example of creating your own collision settings, using the default values:
```julia
collision_settings = CollisionSettings(
  collisions_on = true,
  floe_floe_max_overlap = 0.55,
  floe_domain_max_overlap = 0.75,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.
#### Fracture Settings
`FractureSettings` changes Subzero's floe fractures. The availible settings are:
- fractures_on, which turns this process on and off
- criteria, which sets the rules for which floes fractures (see more below)
- Δt, which sets the number of timesteps between this process running
- deform_on, which turns on a sub-process where a floe is deformed around the largest collision on that timestep before it is fractured
- npieces, which sets the number of pieces a floe should try to fracture into (note that this number might not be met depending on floe shape)
- nhistory, which sets the length of stress history, stress from the previous `nhistory` timesteps, to record in order to determine a floe's current stress, which is the mean of the stress history

Here is an example of creating your own fracture settings, using the default values:
```julia
fracture_settings = FractureSettings(
  fractures_on = false,
  criteria = NoFracture(),
  Δt = 0,
  deform_on = false,
  npieces = 3,
  nhistory = 1000,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

If you do want to turn fractures on, you will need to create a criteria that is not the NoFracture() criteria. Right now, there is only one other fracture criteria, which is Hibler's elliptical yield curve. To learn more about this criteria, see the 1979 paper "A Dynamic Thermodynamic Sea Ice Model." This criteria uses the mean floe height and take in two parameters, `pstar` and `c`, which can be tuned to get the fracture behavior that you want. A simple way to create this criteria is as follows:
```julia
pstar = 2.25e5
c = 20.0
criteria = HiblerYieldCurve(floe_arr, pstar, c)
```
#### Simplification Settings
`SimplificationSettings` changes Subzero's floe simplification. The availible settings are:
- dissolve_on, which turns on and off the dissolving of small floes
- min_floe_area, which sets the minimum size for floes before they are dissolved, and also a barrier for making floes smaller through fracture and other processes
- smooth_vertices_on, which turns on and off the process of smoothing floe vertices to decrease the total number of vertices
- max_vertices, the total number of verticies a floe can have before smoothing
- Δt_smooth, which sets the number of timesteps between floe smoothing

Here is an example of creating your own simplification settings, using the default values:
```julia
simp_settings = SimplificationSettings(
    dissolve_on = true,
    min_floe_area = 1e6,
    smooth_vertices_on = true,
    max_vertices = 30,
    Δt_smooth = 20,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

### Timesteps
You have the ability to set the simulation's timestep in seconds using `∆t` and set the total number of timsteps the simulation will run for, `n∆t`. The default is `∆t = 10` seconds and `n∆t = 7500` timesteps. 

### Output Writers
You can add four types of output writers, and as many of each type as you would like. The four types are as follows: `InitialStateOutputWriter`, `CheckpointOutputWriter`, `FloeOutputWriter`, and `GridOutputWriter`. When any of these objects are create, the file that they will write to is also created automatically. A brief desctiption of each is below:
#### InitialStateOutputWriter
The initial state output writer allows you to save the initial state of your simulation so that you can re-load it later. It saves the simulation object that you are currently creating to a JLD2 file. An `InitialStateOutputWriter` has two fields:
- `filename`, including the path to the file, to save the file to
- `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.

You can create an `InitialStateOutputWriter` directly as a struct or with the following function call:
```julia
init_writer = InitialStateOutputWriter(
    dir = ".",
    filename = "initial_state.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
```
Note that these are the default values and you only need to pass in arguments that you wish to change from these values. Also note that all inputs are keyword arguments.

#### CheckpointOutputWriter
The checkpoint output writer allows you to save the state of the floes, ocean, and atmosphere at a specified part of the simulation into a JLD2 file. This will give you the ability to easily restart your simulation from the last-saved checkpoint if it were to fail for any reason, or if you wanted to continue a previous run from its stopping point. A `CheckpointOutputWriter` has three fields:
- `Δtout`, which specifies the umber of timesteps between checkpoint outputs starting from the first timestep
- `filename`, including the path to the file, to save the file to
-  `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.

You can create an `CheckpointOutputWriter` directly as a struct or with the following function call:
```julia
checkpointer = CheckpointOutputWriter(
    Δtout,
    dir = ".",
    filename = "checkpoint.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
```
Note that other than `∆tout`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

#### FloeOutputWriter
The floe output writer allows you to save floe values into a JLD2 file. This will give you the ability to easily analyze floe fields across timesteps. A `FloeOutputWriter` has four field:
- `outputs`, which specifies which floe fields should be included
- `Δtout`, which specifies the umber of timesteps between floe outputs starting from the first timestep
- `filename`, including the path to the file, to save the file to
- `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.

You can create an `FloeOutputWriter` directly as a struct or with the following function call:
```julia
floewriter = FloeOutputWriter(
    outputs,
    Δtout,
    dir = ".",
    filename = "floes.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
```
The `outputs` field takes in a list of symbols corresponding to floe fields. For example, if you want the floe output writer to output the floes centroid and coordinates then `outputs = [:centroid, :coords]`. If you want all floe fields then you can simply admit the outputs field altogether and all floe fields will be output. Note that other than `outputs` and `∆tout`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

#### GridOutputWriter
The grid output writer allows you to floe values averaged onto a course grid to a NetCDF file. This will give you the ability to easily analyze floe characteristics on a grid. A `GridOutputWriter` has eight field:
- `outputs`, which specifies which floe fields should be included
- `Δtout`, which specifies the umber of timesteps between floe outputs starting from the first timestep
- `filename`, including the path to the file, to save the file to
- `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.
- `xg`, the grid lines in the x-direction of the grid that you would like the calculations done over (doesn't have to be the same as simulation grid)
- `yg`, the grid lines in the y-direction of the grid that you would like the calculations done over (doesn't have to be the same as simulation grid)
- `data`, three-dimensional array that holds grid averaged data prior to writing to file
- `average`, boolean that specifies if the gridded data should be averaged over each timestep between writing to file, or if it should just be calculated prior to writing for that singular timestep (NOT IMPLEMENTED YET).

You can create an `GridOutputWriter` directly as a struct or with the following function call:
```julia
gridwriter = GridOutputWriter(
    outputs,
    Δtout,
    grid,
    dims;
    dir = ".",
    filename = "gridded_data.nc",
    overwrite = false,
    average = false,
)
```
The `outputs` field takes in a list of symbols. To see all possible outputs, call the `get_known_grid_outputs()` function. For example, if you want the grid output writer to output the floe masses and areas averaged on the grid then `outputs = [:mass_grid, :area_grid]`. If you want all possible fields then you can simply admit the outputs field altogether and all grid fields will be output. The `grid` field is the simulation grid, and then `dims` field specifies the dimensions of the grid you would like the output calculated on.

Note that other than `outputs`, `∆tout`, `grid`, and `dims`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

#### OutputWriters
Once you have created all of the types of output writers you need, you must combine them into one `OutputWriters` object that will be a simulation field. 
The `OutputWriters` struct has four fields: `initialwriters`, `floewriters`, `gridwriters`, and `checkpointwriters`. For each you can supply a StructArray of the specified type of writer. If you do not have any writers of a given type, don't provide any. Below I create two example `OutputWriters` objects, assuming that I have already created the following outputwriters: `initwriter1`, `checkpointer1`, `floewriter1`, `floewriter2`, and `gridwriter1`.

```julia
using StructArrays, Subzero

outputwriters1 = OutputWriters(
    initialwriters = StructArray([initwriter1]),
    floewriters = StructArray([floewriter1]),
)

outputwriters2 = OutputWriters(
    initialwriters = StructArray([initwriter1]),
    checkpointwriters = StructArray([checkpointer1]),
    floewriters = StructArray([floewriter1, floewriter2]),
    gridwriters = StructArray([gridwriter1]),
)
```
Here you can see that you can choose which values to supply and that you can supply more than one of each type if desired. You might want to do this is you want different outputs at different timeframes. 

### Reproducibility

### Creating the Simulation
Once you have created all of the above objects, you can combine them to create a `Simulation`. A simulation has quite a few fields, all of which are talked about in more detail above. Since there are so many fields, they are keyword defined, so you must provide a keyword when creating the struct. The only necessary argument is the `model` as everything else has a default value. However, a table of optional elements, and their default values, is as follows:

| Keyword           |    Default Value        | What is it?                                               |
| ----------------- | ------------------------| --------------------------------------------------------- |
| consts            | Constants()             | Physical parameters used in the simulation                |
| rng               | Xoshiro()               | Random number generator - can seed for reproducibility    |
| verbose           | false                   | Flag for printing timesteps and updates during simulation |
| name              | "sim"                   | Name of simulation                                        |
| ∆t                | 10                      | Length of simulation timestep in seconds                  |
| n∆t               | 7500                    | Number of timesteps run in simulation                     |
| coupling_settings | CouplingSettings()      | Settings for coupling during the simulation               |
| collision_settings| CollisionSettings()     | Settings for collisions during the simulation             |
| fracture_settings | FractureSettings()      | Settings for fractures during the simulation              |
| simp_settings     | SimplificationSettings()| Settings for floe simplification during the simulation    |
| writers           | OutputWriters()         | Lists of output writers to be written to during the run   |

Here is an example of how to create a simulation, providing some of these optional inputs given that we have already created the following objects: `my_model`, `my_consts`, `my_coupling_settings`, `my_simp_settings`, and `my_writers`.

```julia
my_simulation = Simulation(
    model = my_model,
    consts = my_consts,
    Δt = 5,  # Timestep of 5 seconds
    coupling_settings = my_coupling_settings,
    simp_settings = my_simp_settings,
    writers = my_writers,
)
```
