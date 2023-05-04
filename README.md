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

Let’s run a basic simulation with stationary floes pushed into a collision boundary by a uniform, zonally flowing ocean. In this simulation, collisions between floes are on by default and we will enable floe fracturing.  

```julia 
grid = RegRectilinearGrid(0, 1e5, 0, 1e5, 1e4, 1e4) 
ocean = Ocean(grid, 0.25, 0.0, 0.0) 
atmos = Atmos(grid, 0.0, 0.0, 0.0) 
domain = Domain( 
  CollisionBoundary(grid, North()), 
  CollisionBoundary(grid, South()), 
  CollisionBoundary(grid, East()),
  CollisionBoundary(grid, West()),
) 
floe_arr = initialize_floe_field(50, [0.7], domain, 0.5, 0.05) 
model = Model(grid, ocean, atmos, domain, floe_arr) 

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

simulation = Simulation( 
  model = model, 
  consts = consts, 
  Δt = Δt, 
  nΔt = 5000, 
  verbose = true, 
  fracture_settings,
)

floewriter = FloeOutputWriter(100, dir = "output/sim", filename = "floes.jld2", overwrite = true)
run!(simulation, [floewriter])
``` 

Check out our documentation for more examples. Below you will find videos from simulations initially described in Manucharyan and Montemuro’s paper, as well as a coupled simulation with Oceananigans.  

## Citing:
TODO

## Contributing:
If you’re interested in helping develop Subzero, have found a bug, or have a new feature that you want implemented, please open an issue on the repository.  For more information, please see our contributor's guide.  

## Movies:

## Performance Benchmarks: 
INSERT RESULTS
