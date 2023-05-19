<!-- Title -->
<h1 align="center">
  Subzero.jl
</h1>

<!-- description -->
<p align="center">
  <strong>:ice_cube: Sea-ice model to explicitly simulate individual floe life cycles using complex discrete elements with time-evolving shapes. Easy to couple with Oceananigans for explorations of coupled ice-ocean dynamics. Ported from MATLAB to Julia. :ice_cube:</strong>
</p>

<!-- Information badges -->
<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  </a>
</p>

<!-- CI/CD badges -->
<p align="center">
  <a href="https://github.com/Caltech-OCTO/Subzero.jl/actions/workflows/CI.yml?query=branch%3Amain">
    <img alt="GitHub Actions CI Status" src="https://github.com/Caltech-OCTO/Subzero.jl/actions/workflows/CI.yml/badge.svg?branch=main">
  </a>
  <a href="https://codecov.io/gh/Caltech-OCTO/Subzero.jl">
    <img alt="CodeCov Status" src="https://codecov.io/gh/Caltech-OCTO/Subzero.jl/branch/main/graph/badge.svg">
  </a>
</p>


Subzero is an easy-to-use Julia translation of Manucharyan and Montemuro’s model described in the paper “SubZero: A Sea Ice Model With an Explicit Representation of the Floe Life Cycle.” The model has been restructured to leverage Julia’s language abstractions for ease of setting up new simulation runs and allowing more types of simulations without code modifications. It has been designed for stand-alone simulations, or to be coupled per-timestep with ocean model Oceananigans.  

Subzero.jl was ported and restructured by Skylar Gering and originally developed by Georgy Manucharyan and Brandon Montemuro.

## Contents

* [Installation instructions](#installation-instructions)
* [Running your first model](#running-your-first-model)
* [Citing](#citing)
* [Contributing](#contributing)
* [Movies](#movies)
* [Performance benchmarks](#performance-benchmarks)

## Installation Instructions:

Subzero is not yet a publically registered Julia package. If you have access to this repository, you have ability to use it as a package added directly from GitHub. This will require having a SSH key on your computer and stored in GitHub so that you are able to securely use this code as it is in private repository for now. 

You can create and add a SSH key using the following instructions from GitHub: 
1. [Checking for existing SSH keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/checking-for-existing-ssh-keys)
2. [Generating a new SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
3. [Adding SSH key to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
4. [Checking your SSH connection](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/testing-your-ssh-connection)

Once you have established your SSH connection from within terminal, you now need to update your Julia `startup.jl` file. This is within the `.julia/config` folder. If you don't have a `.config` folder, please make one using `mkdir config` run on terminal within your `.julia` folder. If you don't have a `startup.jl` file you can make this using `touch startup.jl` within the `config` folder. Then, using a text editor (such as vim), add the following line to your `startup.jl` file: 

`ENV["JULIA_PKG_USE_CLI_GIT"]=true`

We need this as Julia's SSH library can't read the types of SSH keys that GitHub now requires. This will have Julia use your local command line interface (CLI) version of Git. This only works with Julia 1.7 and higher. 

At this point, you can start a Julia REPL and enter into the package manager mode. In package manager mode, you can simply give the command

`add "git@github.com:Caltech-OCTO/Subzero.jl.git"`

which will add Subzero as a package. 

## Running your first model:

For detailed instructions on how to create different types of models and simulations, please see the [documentation.md file](https://github.com/Caltech-OCTO/Subzero.jl/blob/documentation/documentation.md). However, we will give a basic example here.

Let’s run a basic simulation with initially stationary floes pushed into a collision boundary by a uniform, zonally flowing ocean. In this simulation, collisions between floes are on by default and we will enable floe fracturing.  

```julia
const FT = Float64
# Create environment
grid = RegRectilinearGrid(FT, 0, 1e5, 0, 1e5, 1e4, 1e4) 
ocean = Ocean(FT, grid, 0.25, 0.0, 0.0) 
atmos = Atmos(FT, grid, 0.0, 0.0, 0.0) 
# Create domain
domain = Domain( 
  CollisionBoundary(FT, North, grid), 
  CollisionBoundary(FT, South, grid), 
  CollisionBoundary(FT, East, grid),
  CollisionBoundary(FT, West, grid),
)
# Create floes
floe_arr = initialize_floe_field(50, [0.7], domain, 0.5, 0.05) 
# Create model
model = Model(grid, ocean, atmos, domain, floe_arr) 
# Create simulation
modulus = 1.5e3*(mean(sqrt.(floe_arr.area)) + minimum(sqrt.(floe_arr.area))) 
consts = Constants(E = modulus) 
fracture_settings = FractureSettings( 
  fractures_on = true,
  criteria = HiblerYieldCurve(floe_arr),
  Δt = 75,
  npieces = 3,
  nhistory = 1000,
  deform_on = true, 
) 
floewriter = FloeOutputWriter(100, dir = "output/sim", filename = "floes.jld2", overwrite = true)
writers = OutputWriters(floewriters = StructArray([floewriter]))
simulation = Simulation( 
  model = model, 
  consts = consts, 
  Δt = Δt, 
  nΔt = 5000, 
  verbose = true, 
  fracture_settings,
  writers = writers,
)
# Run simulation
run!(simulation)
``` 

Check out our [documentation](https://github.com/Caltech-OCTO/Subzero.jl/blob/documentation/documentation.md) for more examples and explanations for the code above.

## Citing:
We still need to figure this out. Please reach out so we can discuss.

## Contributing:
If you’re interested in helping develop Subzero, have found a bug, or have a new feature that you want implemented, please open an issue on the repository and we can talk about this.  We will be working on a contributers' guide for this in the future.

## Movies:

**Shear Flow**
In this simulation, the ocean flow is 0m/s at the top and bottom of the domain, gradually increasing towards 0.5m/s in the middle of the domain. All four boundaries are periodic. We used a timestep of 20 seconds for 4,320 timesteps, which is one day. 

<img src="https://github.com/Caltech-OCTO/Subzero.jl/assets/60117338/2b13746e-e4db-4ceb-92c5-59f50f2cab32" alt="Shear flow" width="800">

**Simple Strait**
In this simulation, the ocean floe is uniformly -0.3 m/s from top to bottom of the simulation. The top and bottom boundaries are periodic, with the right and left being collision boundaries. However, the collision boundaries are covered by two pieces of topography forming the strait. This simulation also has 20 second timesteps, run for 4,320 timesteps, which is one day.

<img src="https://github.com/Caltech-OCTO/Subzero.jl/assets/60117338/ec331900-aeb7-4b05-a713-a0d5a2b529b8" alt="Simple strait" width="800">

## Performance Benchmarks: 
Here we compare Subzero runtimes in Julia and MATLAB for the shear flow simulation. The code was run on Caltech's HPC cluster. The Julia code was run with 12 threads (12 CPUs per 1 node and 1 task). The MATLAB code has parfor loops with 12 workers (12 tasks, 1 CPU per task).

<img src="https://github.com/Caltech-OCTO/Subzero.jl/assets/60117338/9532e883-0f1d-4399-b713-de24803de72f" alt="Performance data" width="800">

