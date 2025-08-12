# How to Create a Simulation: 

You have a lot of flexibility in designing a simulation within Subzero. First, you need to create all the pieces of your model, including the grid, the ocean, the atmosphere, the domain, and the floes. After than, you can use your model to create a simulation, where you will also specify physical constants and other runtime parameters.

**The documentation, and to some extent the source code, is being cleaned up. This means that right now, some of the documentation is here, and some is in the [tutorial](https://caltech-octo.github.io/Subzero.jl/dev/tutorial/) and [API](https://caltech-octo.github.io/Subzero.jl/dev/api/) sections of the documentation website. Start with the tutorial and API sections and then come back here.**

## Contents
* [Building a Model](#building-a-model)
    * [Grid, Domain, Ocean, Atmosphere](#grid-domain-ocean-atmosphere)
    * [Floes](#floes)
    * [Making the Model](#making-the-model)
* [Building a Simulation](#building-a-simulation)
    * [Constants](#constants)
    * [Physical Process Settings](#physical-process-settings)
    * [Timesteps](#timesteps)
    * [Output Writers](#output-writers)
    * [Reproducibility](#reproducibility)
    * [Creating the Simulation](#creating-the-simulation)

## Building a Model

See the [tutorial](https://caltech-octo.github.io/Subzero.jl/dev/tutorial/) and [API](https://caltech-octo.github.io/Subzero.jl/dev/api/) sections of the documentation website for information on how to setup a `Model`.

Start with these and then return here to learn how to make the rest of a simulation!

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

In particular, Young's Modulus is usually calculated using the total floe area after floe initialization in the original Subzero code:
```julia
E = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area)))
```

### Physical Process Settings
Subzero allows you to turn on and off various physical processes, as well as change the settings for these physical processes. The physical processes availible are: coupling, collisions, fractures, and simplfication. We will be adding corner fracturing, packing, rafting, ridging, and welding. For each of these physical processes, you can create a settings object, and change the parameters within that object to change the physical process to meet your simulation needs. 
#### Coupling Settings
`CouplingSettings` changes how Subzero uses the ocean and atmosphere two-dimensional vector fields within the `Ocean` and `Atmos` structs, whether those values are static user-provided values or values that are updated every timesteps through coupling with Oceananigans. The availible settings are:
- coupling_on, which turns this process on and off
- ∆t, which sets the number of timesteps between this process running
- ∆d, which sets the number of ocean/atmosphere grid cells around a floe to consider when interpolating ocean/atmosphere forcings
- two_way_coupling_on, which turns on and off the calculation of the ice and atmosphere effects on the ocean and stores output in the ocean's stress fields.

Here is an example of creating your own coupling settings, using the default values:
```julia
couple_settings = CouplingSettings(
  coupling_on = true,
  Δt = 10,
  Δd = 1,
  two_way_coupling_on = false,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

#### Collision Settings
`CollisionSettings` changes Subzero's floe-floe and floe-domain collisions. The availible settings are:
- collisions_on, which turns this process on and off
- floe_floe_max_overlap, which sets the maximum fraction of floe-floe overlap before fusing the two floes together 
- floe_domain_max_overlap, which sets the maximum fraction for floe-domain overlap before removing the floe from the simulation

Here is an example of creating your own collision settings, using the default values:
```julia
collision_settings = CollisionSettings(
  FT,
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
- deform_on, which turns on a sub-process where a floe is deformed around the floe that is overlaps with the most during that timestep before it is fractured
- npieces, which sets the number of pieces a floe should try to fracture into (note that this number might not be met depending on floe shape)

Here is an example of creating your own fracture settings, using the default values:
```julia
fracture_settings = FractureSettings(
  fractures_on = false,
  criteria = NoFracture(),
  Δt = 0,
  deform_on = false,
  npieces = 3,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

If you do want to turn fractures on, you will need to create a criteria that is not the `NoFracture()` criteria. Right now, there is are two other fracture criteria.

One is Hibler's elliptical yield curve. To learn more about this criteria, see the 1979 paper "A Dynamic Thermodynamic Sea Ice Model." This criteria uses the mean floe height and take in two parameters, `pstar` and `c`, which can be tuned to get the fracture behavior that you want. The `HilberYieldCurve` object then has both of these parameters as fields as well as the `vertices` of the yield curve. A simple way to create this criteria is as follows:
```julia
pstar = 2.25e5
c = 20.0
criteria = HiblerYieldCurve(FT, floe_arr, pstar, c)
```

The other is Mohr's Cone yield curve. To learn more about this critera, see the 2009 paper "Coulombic faulting from the grain scale to the geophysical scale: lessons from ice." This criteria takes in three parameters that define the shape of the cone. They are as follows:
- q: based on the coefficient of internal friction (µi) by (μi^2 + 1)^(1/2) + μi^2
- σc: uniaxial compressive strength
- σ11: negative of the x-coordinate of one vertex of cone (triangle in 2D) and negative of the y-coordinate of adjacend vertex in principal stress space

Note that the yield curve does not depend on the floe field. The `MohrsCone` object then has a `vertices` field that holds the coordiantes of the cone in principle stress space.


#### Simplification Settings
`SimplificationSettings` changes Subzero's floe simplification. The availible settings are:
- smooth_vertices_on, which turns on and off the process of smoothing floe vertices to decrease the total number of vertices
- max_vertices, the total number of verticies a floe can have before smoothing
- Δt_smooth, which sets the number of timesteps between floe smoothing
- tol, which is the tolerance in Douglas–Peucker's polygon simplification algorithm in meters.

Here is an example of creating your own simplification settings, using the default values:
```julia
simp_settings = SimplificationSettings(
    FT,
    min_floe_area = 1e6,
    smooth_vertices_on = true,
    max_vertices = 30,
    Δt_smooth = 20,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

#### Ridging and Rafting Settings
`RidgeRaftSettings` changes Subzero's ridging and rafting process. The availible settings are:
- ridge_raft_on, which turns on and off the process of smoothing floe vertices to decrease the total number of vertices
- Δt  sets the number of timesteps between floe smoothing
- ridge_probability is the likelyhood (between 0-1) that two floes that meet the criterion actually ridge
- raft_probability is the likelyhood (between 0-1) that two floes that meet the criterion actually raft
- min_overlap_frac is the minimum overlap area fraction between a floe and another floe/domain for that floe to ridge or raft
- min_ridge_height is the minimum floe height to ridge with a floe/domain
- max_floe_ridge_height is the maximum floe height to ridge with another floe
- max_domain_ridge_height is the maximum floe height to ridge with a domain element
- max_floe_raft_height is the maximum floe height to raft with another floe
- max_domain_raft_height is the maximum floe height to raft with a domain element
- domain_gain_probability is the probalility that a floe that rafts with a domain element keeps all of its mass (0) or if that mass is removed and lost to the domain element (1).

Here is an example of creating your own simplification settings, using the default values:

```julia
ridgeraft_settings = RidgeRaftSettings(
  ridge_raft_on = false,
  Δt = 0,
  ridge_probability = 0.95,
  raft_probability = 0.95,
  min_overlap_frac = 0.01,
  min_ridge_height = 0.2,
  max_floe_ridge_height = 5.0,
  max_domain_ridge_height = 1.25,
  max_floe_raft_height = 0.25,
  max_domain_raft_height = 0.25,
  domain_gain_probability = 1.0,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

#### Welding Settings
`WeldSettings` changes Subzero's welding process. The availible settings are:
- weld_on is a boolean flag for if welding should be turned on in the simulation
- Δts is a a list of multiples of timesteps during which welding code will run, welding will be run at multiples of all elements, each with domain split into corresponding Nx and Ny values
- Nxs is a list of number of x-directional bins to split the domain into at corresponding timesteps
- Nys is a list of number of x-directional bins to split the domain into at corresponding timesteps
- min_weld_area is the minimum area a weld can create for two floes to weld
- max_weld_area is the maximum area a weld can create for two floes to weld
- welding_coeff is a non-dimensional parameter, multiplied by ratio of overlap between two floes to original floe area to determin probability that a floe will merge. The larger this is, the more likely floes are to weld. Probability with 5% overlap is `welding_coeff * (0.05) > rand()`

Here is an example of creating your own welding settings, using the default values:

```julia
weld_settings = WeldSettings(
    weld_on = false,
    Δts = Vector{Int}(),
    Nxs = Vector{Int}(),
    Nys = Vector{Int}(),
    min_weld_area = 1e6,
    max_weld_area = 2e9,
    welding_coeff = 150,
)
```
Note that you only need to provide values to the fields that you wish to change from the defaults.

### Timesteps
You have the ability to set the simulation's timestep in seconds using `∆t` and set the total number of timsteps the simulation will run for, `n∆t`. The default is `∆t = 10` seconds and `n∆t = 7500` timesteps. 

### Output Writers
You can add four types of output writers, and as many of each type as you would like. The four types are as follows: `InitialStateOutputWriter`, `CheckpointOutputWriter`, `FloeOutputWriter`, and `GridOutputWriter`. When any of these objects are created, the file that they will write to is also created automatically. A brief desctiption of each is below:
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
- `Δtout`, which specifies the umber of timesteps between floe outputs starting from the first timestep
- `outputs`, which specifies which floe fields should be included
- `filename`, including the path to the file, to save the file to
- `overwrite` boolean that specifies whether a file with the same filename should be overwritter or if an error should be thrown during creation.

You can create an `FloeOutputWriter` directly as a struct or with the following function call:
```julia
floewriter = FloeOutputWriter(
    Δtout;
    outputs = collect(fieldnames(Floe)),
    dir = ".",
    filename = "floes.jld2",
    overwrite = false,
    jld2_kw = Dict{Symbol, Any}(),
)
```
The `outputs` field takes in a list of symbols corresponding to floe fields. For example, if you want the floe output writer to output the floes centroid and coordinates then `outputs = [:centroid, :coords]`. If you want all floe fields then you can simply omit the outputs field all together and all floe fields will be output. Note that other than `∆tout`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

Note that if you have Periodic calls, and thus ghost floes in your simulation, these will also be saved by the `FloeOutputWriter`. If you want to exclude these floes from your analysis or when otherwise using the `FloeOutputWriter` output, you can do so by only including floes with a `ghost_id = 0`.

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
    FT,
    Δtout,
    grid,
    dims;
    outputs = collect(get_known_grid_outputs()),
    dir = ".",
    filename = "gridded_data.nc",
    overwrite = false,
    average = false,
)
```
The `outputs` field takes in a list of symbols. To see all possible outputs, call the `get_known_grid_outputs()` function. For example, if you want the grid output writer to output the floe masses and areas averaged on the grid then `outputs = [:mass_grid, :area_grid]`. If you want all possible fields then you can simply admit the outputs field altogether and all grid fields will be output. The `grid` field is the simulation grid, and then `dims` field specifies the dimensions of the grid you would like the output calculated on.

Note that other than `∆tout`, `grid`, and `dims`, all values have default values so you only need to pass in arguments that you wish to change from these values. Furthermore, all inputs but ∆tout are keyword arguments, so you must use the keyword when passing in new values.

#### OutputWriters
Once you have created all of the types of output writers you need, you must combine them into one `OutputWriters` object that will be a simulation field. 
The `OutputWriters` struct has four fields: `initialwriters`, `floewriters`, `gridwriters`, and `checkpointwriters`. For each you can supply a StructArray of the specified type of writer. If you do not have any writers of a given type, don't provide any. Below I create two example `OutputWriters` objects, assuming that I have already created the following outputwriters: `initwriter1`, `checkpointer1`, `floewriter1`, `floewriter2`, and `gridwriter1`.

```julia
using StructArrays, Subzero

outputwriters1 = OutputWriters(initwriter1, floewriter1)

outputwriters2 = OutputWriters(
   initwriter1,
   checkpointer1,
   floewriter1,
   floewriter2,
   gridwriter1,
)
```
Here you can see that you can choose which values to supply and that you can supply more than one of each type if desired. You might want to do this is you want different outputs at different timeframes. 

### Reproducibility
The simulations are currently completly reproducible when run single-threaded. The simulation takes a `rng` argument, and if it is provided with a seeded random number generator, the two simulations with the same set of starting floes will produce the exact same results. Note that the same set of floes can be reproduced by providing a seeded random number generator to the floe creation functions. 

This has not yet been achieved for multi-threaded runs. 

### Creating the Simulation
Once you have created all of the above objects, you can combine them to create a `Simulation`. A simulation has quite a few fields, all of which are talked about above in more detail. Since there are so many fields, they are keyword defined, so you must provide a keyword when creating the struct. The only necessary argument is the `model` as everything else has a default value. However, a table of optional elements, and their default values, is as follows:

| Keyword           |    Default Value        | What is it?                                               |
| ----------------- | ------------------------| --------------------------------------------------------- |
| consts            | Constants()             | Physical parameters used in the simulation                |
| rng               | Xoshiro()               | Random number generator - can seed for reproducibility    |
| verbose           | false                   | Flag for printing timesteps and updates during simulation |
| name              | "sim"                   | Name of simulation                                        |
| ∆t                | 10                      | Length of simulation timestep in seconds                  |
| n∆t               | 7500                    | Number of timesteps run in simulation                     |
| floe_settings     | FloeSettings()          | Settings for making new floes    during the simulation    |
| coupling_settings | CouplingSettings()      | Settings for coupling during the simulation               |
| collision_settings| CollisionSettings()     | Settings for collisions during the simulation             |
| fracture_settings | FractureSettings()      | Settings for fractures during the simulation              |
| simp_settings     | SimplificationSettings()| Settings for floe simplification during the simulation    |
| ridgeraft_settings| RidgeRaftSettings()     | Settings for ridging and rafting during the simulation    |
| weld_settings     | WeldSettings()          | Settings for welding during the simulation                |
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

### Running the simulation

You can now use the `run!` function to run the simulation:
```
run!(my_simulation)
```

If you wish to couple to Oceananigans, you will need to run each model timestep by timestep and pass the needed fields back and forth. You can run a single timestep of the simulation using the `timestep_sim!` function. This also needs the current timestep the simulation is on (`tstep`) as an argument.

```
timestep_sim!(
   my_simulation,
   tstep,
)
```

Note that we are working on a more elegant solution to coupling with Oceananigans and CliMA and this page will be updated once that is in place. 

If you run your simulation in multiple parts and need to re-start your simulation from files, the `restart!` function will be a good place to start. However, note that it is quite simple and users may need to write their own restart function if they want any complex behavior. 

The provided `restart!` function takes in the output file from both an `InitialStateOutputWriter` and a `CheckpointOutputWriter` to restart the simulation. In addition to providing these two files, the user must also provide the number of timesteps to run the next part of the simulation for (`new_nΔt`) and new output writers. The user also has an option to specify a non-zero starting timestep for the simulation using the keyword argument `start_tstep`.

```
restart!(
   initial_state_fn,
   checkpointer_fn,
   new_nΔt,
   new_output_writers;
   start_tstep = 0,
)
```

### Plotting

If your simulation has both a `FloeOutputWriter` and an `InitialStateOutputWriter`, you can use the built in plotting function to make an MP4 file with each frame as a timestep saved by the `FloeOutputWriter`. You do this as follows:

```
plot_sim(
    floe_output_writer_file_path,
    initial_state_output_writer_file_path,
    Δt,
    output_file_path,
)
```

where `floe_output_writer_file_path` is the .jl file saved by the `FloeOutputWriter`, `initial_state_output_writer_file_path` is the .jl file saved by the `InitialStateOutputWriter`, and `output_file_path` is the file path and name you want your .mp4 file saved as. `Δt` is the model timestep.
