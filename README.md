<!-- Title -->
<h1 align="center">
  Subzero.jl
</h1>

<!-- description -->
<p align="center">
  <strong>:ice_cube: UW's Sea-ice model SubZero ported from MATLAB to Julia :ice_cube:</strong>
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


As of now, Subzero is not a registered Julia package for privacy reasons. However, it can still be added like a package. To do this, you will need to add a SSH key on your computer. Instructions for that can be found [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent). When you follow step 2, add `-m PEM`, as in `$ ssh-keygen -m PEM -t ed25519 -C "your_email@example.com"`. Then continue on to [add the key to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account) and [check your connection](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/testing-your-ssh-connection).

After this, if you are a member of this repository, you can add Subzero as a package by writing this in the package manager: 

`add git@github.com:Caltech-OCTO/Subzero.jl.git`. 

If you want to work with a specific branch of the repository use: 

`add "git@github.com:Caltech-OCTO/Subzero.jl.git"#branch_name`. 

Note the quotes in this second option.

Once you you have added Subzero as a package, you now just need to put `using Subzero` at the top of any scripts using Subzero. 

For example scripts, see the examples folder. 
