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

Subzero is not yet a registered Julia package. If you have access to this repository, you have ability to use it as a package from our private package registry.  NEED TO SETUP 

## Running your first model:

Let’s run a basic simulation with stationary floes pushed into a collision boundary by a uniform, zonally flowing ocean. In this simulation, collisions between floes are on by default and we will enable floe fracturing.  

```julia 
grid = RegRectilinearGrid(Float64, 0, 1e5, 0, 1e5, 1e4, 1e4) 
ocean = Ocean(Float64, grid, 0.25, 0.0, 0.0) 
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
