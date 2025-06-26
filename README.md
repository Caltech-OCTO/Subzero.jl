![Subzero.jl](https://github.com/Caltech-OCTO/Subzero.jl/blob/main/docs/src/assets/title.gif?raw=true)

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://caltech-octo.github.io/Subzero.jl/dev/)
[![CI](https://github.com/Caltech-OCTO/Subzero.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Caltech-OCTO/Subzero.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![CodeCov](https://codecov.io/gh/Caltech-OCTO/Subzero.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Caltech-OCTO/Subzero.jl)
[![Status](https://www.repostatus.org/badges/latest/active.svg?style=flat-square)](https://www.repostatus.org/#active)


## _Fast and Flexible Sea Ice Dynamics_

Subzero.jl is a native [Julia](https://julialang.org/) discrete-element model (DEM) for exploring fine-scale sea ice dynamics, reimplementing and enhancing MATLAB model [SubZero by Manucharyan and Montemuro](https://doi.org/10.1029/2022MS003247).

- 🚀 Runs over **35 times faster** that original MATLAB model for title simulation!
- 🧩 Modular simulation model makes it easy to **customize simulations**!
  - Enable and disable physical processes such as fracturing, ridging, and welding
  - Choose algorithms for key processes (or add your own!)

## _Documentation_

To learn how to build and run simulations, [check out our documentation and tutorials](https://caltech-octo.github.io/Subzero.jl/dev/)!

## _Installation_

Subzero is not yet a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/) (version 1.9 or later). We recommend using [`JuliaUp`](https://github.com/JuliaLang/juliaup) so it is easy to change versions in the future.

2. As `Subzero.jl` is not yet registered as an official package, you will need to install it from GitHub. To do this, you will need a SSH key on your computer and stored in GitHub. GitHub provides documentation for the needed steps: 
    - [Checking for existing SSH keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/checking-for-existing-ssh-keys)
    - [Generating a new key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
    - [Adding a new key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
    - [Testing your connection](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/testing-your-ssh-connection).

3. Once you have established your SSH connection from within terminal, you now need to update your Julia `startup.jl` file. This is within the `.julia/config` folder. If you don't have a `.config` folder, please make one using `mkdir config` run on terminal within your `.julia` folder. If you don't have a `startup.jl` file you can make this using `touch startup.jl` within the `config` folder. Then, using a text editor (such as vim), add the following line to your `startup.jl` file:

    `ENV["JULIA_PKG_USE_CLI_GIT"]=true`

    We need this as Julia's SSH library can't read the types of SSH keys that GitHub now requires. This will have Julia use your local command line interface (CLI) version of Git. This only works with Julia 1.7 and higher. 

4. Launch Julia and enter into the package manager mode by typing `]` in the terminal.
   
5. Run the following

    ```julia
    pkg> add "git@github.com:Caltech-OCTO/Subzero.jl.git"
    ```

    This will add the package to your package manager. After that you return to the REPL mode by hitting the backspace and you are ready to use Subzero! 

    ```julia
    julia> using Subzero
    ```

## _Contributing_

If you’re interested in contributing to the development Subzero, we would love to have you! We welcome all kinds of contributions from bug reports, to documentation, to features, and suggestions. **We can't wait to talk to you.**

Please see [CONTRIBUTING](CONTRIBUTING.md) for more details.

## _Citing_

If you use Subzero.jl as part of your research, teaching, or other activities, we would be grateful if you could cite our work.
We are currently working on a JOSS paper, which will be linked here. If you are ready to publish before that, please reach out to us to discuss citations.

## _Authors_

- Primary Author: [**Skylar Gering (@skygering)**](https://github.com/skygering)

The list of [Subzero contributors](https://github.com/Caltech-OCTO/Subzero.jl/graphs/contributors):

[![Contributers](https://contrib.rocks/image?repo=Caltech-OCTO/Subzero.jl)](https://github.com/Caltech-OCTO/Subzero.jl/graphs/contributors)
