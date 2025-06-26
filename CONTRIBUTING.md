# Subzero.jl contributor guide

Welcome to the Subzero.jl contributor guide. Here you can find information on:

- [Local Repository](#local-repository)
- [Documentation](#documentation)
- [Reporting issues](#reporting-issues)
- [Code changes](#code-changes)

If you are new to open source development, here is a good guide to get started: [first-contributions][first-contributions]. If you want something Julia specific, check out this video: [Open source, Julia
packages, git, and GitHub][tim-git].

## Local Repository
If you want to work on developing Subzero (other than for small documentation changes), you probably want to work locally on your computer. For this, you will want to create a fork of the repository, and then eventually open a pull request to the original code base. Here is a good step-by-step guide of this process: [collaborate on an existing package][existing-package].

## Documentation

Contributing to the documentation is a great way to get involved in Subzero development. If something in the documents is confusing, please feel free to fix it! We always welcome improved documentation. 

Small changes can be done easily in GitHub's web interface (see [Editing
files][gh-edit-files]). Every page in the documentation have an `Edit on GitHub` button at
the top, which takes you to the correct source file. The video [Making Julia documentation
better][tim-doc] guides you through these steps.

If you want to edit larger sections of the documentation, you probably want to use the above instructions to work on a [local repository](#local-repository). You can then make changes there and open a pull request to have them incorporated into the package. 

To see your changes live before opening a pull request, you can locally build the documentation. For this, you will need two terminal windows. Start by making sure you are in the `Subzero/docs` folder in both terminals. Then, in both terminals, you will want to launch Julia.

In the first terminal:
```julia
pkg> activate .
julia> include("make.jl")
```

This will build the documents. It might take a little while.

```julia
pkg> activate .
julia> using LiveServer
julia> serve(;dir = "build")
```

This uses [`LiveServer.jl`][liveserver] to launch a local webserver which you can visit at
[http://localhost:8000](http://localhost:8000). 

For more information on documentation see:

**Useful resources**
 - General information about documenting Julia code in the [Julia manual][julia-doc].
 - [Documentation for `Documenter.jl`][documenter] which is used to render the HTML pages.
 - [Documentation for `Literate.jl`][literate] which is used for tutorials/examples.

## Reporting issues


[first-contributions]: https://github.com/firstcontributions/first-contributions
[tim-git]: https://youtu.be/cquJ9kPkwR8
[gh-edit-files]: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files#editing-files-in-another-users-repository
[tim-doc]: https://youtu.be/ZpH1ry8qqfw
[existing-package]: https://www.matecdev.com/posts/julia-package-collaboration.html
[liveserver]: https://github.com/tlienart/LiveServer.jl
[julia-doc]: https://docs.julialang.org/en/v1/manual/documentation/
[documenter]: https://juliadocs.github.io/Documenter.jl/
[literate]: https://fredrikekre.github.io/Literate.jl/v2/
