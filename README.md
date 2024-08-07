<!-- Title -->
<h1 align="center">
  Subzero.jl
</h1>

<!-- description -->
<p align="center">
  <strong> Fast and Flexible Sea Ice Dynamics </strong>
</p>

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

Subzero.jl is a native [Julia](https://julialang.org/) discrete-element sea ice model for exploring fine-scale sea ice dynamics, 
reimplementing MATLAB model [SubZero by Manucharyan and Montemuro](https://doi.org/10.1029/2022MS003247).

Subzero.jl is **fast** and **flexible**.

- Runs over **35 times faster** that original MATLAB model for title shear flow simulation!
- Modular simulation model makes it easy to customize simulations within a _single_ run file!
  - Enable and disable various physical processes such as fracturing, ridging, and welding.
  - Choose from various algorithms for key physical processes (or add your own!)

<br>

## Installation instructions

Subzero is a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/) (version 1.9 or later).

2. Launch Julia and type

```julia
julia> using Pkg

julia> Pkg.add("Subzero")
```

This installs the latest version that's _compatible with your current environment_.

**You're now ready to use Subzero!**

To learn how to build and run simulations, [check out our documentation and tutorials](https://caltech-octo.github.io/Subzero.jl/dev/)!

<br>

## Citing

If you use Subzero.jl as part of your research, teaching, or other activities, we would be grateful if you could cite our work.
We are currently working on a JOSS paper, which will be linked here. If you are ready to publish before that, please reach out to us to discuss citations.

<br>

## Contributing

If you’re interested in contributing to the development Subzero, we would love to have you!
If you have found a bug, or have a new feature that you want to help implemented, please open an issue on the repository. We can't wait to talk to you.

Please also check out [our contributers' guide](https://caltech-octo.github.io/Subzero.jl/dev/contribute/).

<br>

## ✏️ Authors

- Primary Author: [**Skylar Gering (@skygering)**](https://github.com/skygering)

See also the list of [contributors](https://github.com/Caltech-OCTO/Subzero.jl/graphs/contributors) to Subzero.

<a href="https://github.com/Caltech-OCTO/Subzero.jl/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=Caltech-OCTO/Subzero.jl" />
</a>

Made with [contrib.rocks](https://contrib.rocks).
