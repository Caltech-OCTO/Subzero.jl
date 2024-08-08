<!-- Title -->
![Subzero.jl](https://github.com/Caltech-OCTO/Subzero.jl/blob/main/docs/src/assets/title.gif)

<!-- Repo badges -->
<p align="center">
  <a href="https://caltech-octo.github.io/Subzero.jl/dev/">
    <img alt="Docs latest" src="https://img.shields.io/badge/docs-latest-blue.svg" />
  </a>  
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  </a>
  <a href="https://github.com/Caltech-OCTO/Subzero.jl/actions/workflows/CI.yml?query=branch%3Amain">
    <img alt="GitHub Actions CI Status" src="https://github.com/Caltech-OCTO/Subzero.jl/actions/workflows/CI.yml/badge.svg?branch=main">
  </a>
  <a href="https://codecov.io/gh/Caltech-OCTO/Subzero.jl">
    <img alt="CodeCov Status" src="https://codecov.io/gh/Caltech-OCTO/Subzero.jl/branch/main/graph/badge.svg">
  </a>
</p>

<!-- description -->

## Fast and Flexible Sea Ice Dynamics

Subzero.jl is a native [Julia](https://julialang.org/) discrete-element model (DEM) for exploring fine-scale sea ice dynamics, 
reimplementing MATLAB model [SubZero by Manucharyan and Montemuro](https://doi.org/10.1029/2022MS003247).

- 🚀 Runs over **35 times faster** that original MATLAB model for title simulation!
- 🧩 Modular simulation model makes it easy to **customize simulations**!
  - Enable and disable physical processes such as fracturing, ridging, and welding
  - Choose algorithms for key processes (or add your own!)

## Documentation

To learn how to build and run simulations, [check out our documentation and tutorials](https://caltech-octo.github.io/Subzero.jl/dev/)!

## Installation

Subzero is a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/) (version 1.9 or later).

2. Launch Julia and type

```julia
julia> using Pkg

julia> Pkg.add("Subzero")
```

## Citing

If you use Subzero.jl as part of your research, teaching, or other activities, we would be grateful if you could cite our work.
We are currently working on a JOSS paper, which will be linked here. If you are ready to publish before that, please reach out to us to discuss citations.

## Contributing

If you’re interested in contributing to the development Subzero, we would love to have you!
If you have found a bug or have a new feature that you want to help implement, please open an issue on the repository. **We can't wait to talk to you.**

Please also check out [our contributers' guide](https://caltech-octo.github.io/Subzero.jl/dev/contribute/).

## Authors

- Primary Author: [**Skylar Gering (@skygering)**](https://github.com/skygering)

The list of [Subzero contributors](https://github.com/Caltech-OCTO/Subzero.jl/graphs/contributors):

<a href="https://github.com/Caltech-OCTO/Subzero.jl/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=Caltech-OCTO/Subzero.jl" />
</a>
