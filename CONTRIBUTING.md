# Subzero.jl contributor guide

Welcome to the Subzero.jl contributor guide. Here you can find information on:

- [Local Repository](#local-repository)
- [Documentation](#documentation)
- [Reporting issues](#reporting-issues)
- [Code changes](#code-changes)
- [Acknowledgments](#acknowledgments)

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

If you have found a bug or a problem with Subzero you can open an [issue][new-issue]. Try
to include as much information about the problem as possible and some code that
can be copy-pasted to reproduce it (see [How to create a Minimal, Reproducible
Example][so-mre]).

If you can identify a fix for the bug you can submit a pull request without first opening an
issue, see [Code changes](#code-changes).

## Code changes

Bug fixes and improvements to the code, or to the unit tests are always welcome. If you have
ideas about new features or functionality it might be good to first open an
[issue][new-issue] to get feedback before spending too much time implementing something.

When you are ready to make changes, check out the developer docs section of the documentation for insight into how the code is written and organized.

Remember to always include (when applicable): unit tests which exercises the new code,
and updated documentation.

## Acknowledgments

Thank you to the folks at [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) for having such amazing contributor docs. I took lots of inspiration and links from their wonderful, comprehensive page! Check them out for examples of documentation and guides done right!


[first-contributions]: https://github.com/firstcontributions/first-contributions
[tim-git]: https://youtu.be/cquJ9kPkwR8
[gh-edit-files]: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files#editing-files-in-another-users-repository
[tim-doc]: https://youtu.be/ZpH1ry8qqfw
[existing-package]: https://www.matecdev.com/posts/julia-package-collaboration.html
[liveserver]: https://github.com/tlienart/LiveServer.jl
[julia-doc]: https://docs.julialang.org/en/v1/manual/documentation/
[documenter]: https://juliadocs.github.io/Documenter.jl/
[literate]: https://fredrikekre.github.io/Literate.jl/v2/
[new-issue]: https://github.com/Caltech-OCTO/Subzero.jl/issues/new
[so-mre]: https://stackoverflow.com/help/minimal-reproducible-example
