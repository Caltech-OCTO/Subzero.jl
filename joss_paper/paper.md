---
title: 'Subzero.jl: A Coupled, Julia Version of a Discrete-Element Sea Ice Model'
tags:
  - Julia
  - oceanography
  - sea ice
  - dynamics
authors:
  - name: Skylar Gering
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: "1, 2"
  - name: Mukund Gupta
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Samuel Brenner
    orcid: 0000-0002-0826-1294
    affiliation: 1
  - name: Andrew Thompson
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: California Institute of Technology, USA
   index: 1
 - name: Massachusetts Institute of Technology, USA
   index: 2
 - name: Delft Institute of Technology, NL
   index: 3
date: 1 January 2026
bibliography: paper.bib
---

# Summary

Subzero.jl is a discrete-element model (DEM) for simulating the dynamics of two-dimensional polygonal sea-ice floes. The model is a port of the MATLAB-based SubZero [@Montemuro2023] to the Julia programming language [@Bezanson:2017], with substantial re-engineering to increase performance, improve extensibility and usability, and provide native infrastructure for two-way coupling with the external ocean dynamics model Oceananigans.jl [@Ramadhan2020,@Wagner2025]. Subzero.jl is available as a registered Julia package.

# Statement of need

Arctic sea ice extent and concentration continue to decline at rates that are commonly underestimated by climate projection models [@cite]. Potential sources of uncertainty arise from an inaccurate representation of interactions between the ocean and sea ice within these climate models, as well as the assumption of continuum sea-ice dynamics for the sake of reducing computational complexity.

Discrete-element models (DEMs), where each piece of sea ice is represented as an individual simulation element, all of which can dynamically interact, are used to explore these uncertainties and study fine-scale sea ice dynamics. While a range of sea ice DEMs are available [e.g., @Hopkins2004,@Herman2013,@Rabatel2015,@Damsgaard2018, and others], many of these models make a range of simplifications, particularly in the geometric representation of the floes.

The MATLAB-based SubZero model [@Montemuro2023] allows for complex, and evolving, floe shapes but is computationally expensive. Moreover, due to the highly connected nature of sea ice floes and ocean processes [@Horvat2018,@Gupta2022,@Gupta2024,@Brenner2023c], there is a need for DEM simulations to be coupled to a dynamic ocean model to explore two-way feedbacks, which is a feature not readily available for most extant models. 

We present Subzero.jl, a native Julia [@Bezanson2017] version of the MATLAB discrete-element model SubZero [@Montemuro2023, @Manucharyan2022] that addresses both of these two problems. 
# Functionality

Subzero.jl, represents sea ice floes as polygonal elements that move in response to forcing from the atmosphere and the ocean [@Manucharyan2022]. Over time, a given floe may change its horizontal shape and vertical thickness, and fracture into multiple pieces, as a result of interactions with other floes and topographical elements.  SubZero.jl improves upon the MATLAB version in three major ways: (i) a modular interface that limits the need for users to modify source code,, (ii) the ability to couple to a Julia-based ocean model, and (iii) enhancement in computational speed.

## Modular interface

This new Julia implementation allows a more extendable, user-friendly interface relative to the previous MATLAB version. SubZero.jl is designed using various object ‘types’ that are embedded in a tree-structure. The topmost parent is the ‘Simulation’ object, which contains the ‘Floes’, ‘Domain’, ‘Atmospheric forcing’, ‘Ocean forcing’, ‘Physical settings’ and ‘Output writer’ structs. Using dedicated run scripts, the user may modify parameters from these objects and easily turn on/off physical processes such as fracturing, ridging, and welding. Moreover, the use of Julia’s multiple dispatch allows flexibility to add more functionality without changing the overlaying model structure. This flexibility is highlighted in the tutorial and examples provided with the documentation. Additionally, Subzero.jl includes a suite of unit and integration tests that verify consistency with the original MATLAB model and ensure conservation of energy and momentum. These changes allow for more dynamic and flexible usage and development than in the MATLAB version. 
## Ice-Ocean Coupling

The new code now includes a framework that enables two-way coupled simulations with Julia ocean model Oceananigans.jl [@Ramadhan2020,@Wagner2025]. In this configuration, sea-ice floes in Subzero.jl are forced by gridded ocean velocity fields supplied by Oceananigans.jl and interpolated on a mesh carried by each floe. In turn, floe-resolved ice-ocean stresses computed by Subzero.jl are interpolated on the ocean grid and returned to the ocean model as a spatially heterogeneous surface and temporally evolving surface boundary condition. The exchange of fields is implemented using a callback functionality in Oceananigans.jl, allowing serial coupling at user-defined frequencies. This coupled DEM-LES system resolves the mechanical interaction between evolving ocean fields and individual sea ice floes, enabling investigation of floe-scale modulation of upper-ocean dynamics [e.g., @Brenner2025]. An example coupled simulation is demonstrated in  \autoref{fig:oceananigans_coupled_example}.


\begin{figure}
  \includegraphics[width=\linewidth]{oceananigans_coupled_example.png}
  \caption{Example Subzero-Oceananigans coupled simulation. Left: surface fields of ocean vorticity, overlain by sea ice floes, each coloured by their vorticity; Right: the same ocean vorticity field with floes not plotted–the impacts of ice-ocean coupling are still evident in the “patchiness” of the field.}
  \label{fig:oceananigans_coupled_example}
\end{figure}

## Performance improvements
Subzero.jl gives significant performance improvements relative to the original MATLAB implementation of the model.
We used the `shear_flow` example from the repository ([link](https://caltech-octo.github.io/Subzero.jl/dev/examples/shear_flow/)) as a basis for testing speed enhancements.
The simulation includes 50 floes at 75% sea ice concentration, forced by sheared ocean currents in a doubly-periodic domain; we also created a version with 1000 floes, and then implemented both versions in the MATLAB SubZero codebase.
We tested the models on the California Institute of Technology's "Resnick" High Performance Computing system, on a single compute node using 1–16 CPU cores, with each configuration repeated five times.
Results in \autoref{fig:speed_comparison} show comparisons of both end-to-end runtime (including initialization and I/O) and simulation time (time spent advancing the model state).  


\begin{figure}
  \includegraphics[width=\linewidth]{speed_comparison.png}
  \caption{Add caption...}
\end{figure}


Subzero.jl achieves substantial speedups relative to the MATLAB implementation: 7.2–13.0 times faster (end-to-end) and 17.4–31.1 times faster (simulation) for 50 floes, increasing to 32.2–48.5× (end-to-end) and 43.9–61.1 times (simulation) for 1000 floes.
Parallel scaling is modest in both implementations, with simulation time speedups of 1.1 times (Julia) and 1.8 times (MATLAB) for 50 floes, and 1.9 times (Julia) and 2.3 times (MATLAB) for 1000 floes when increasing from 1 to 16 cores. 

The speed increase of Subzero.jl relative to the original MATLAB model will enable longer and more complex simulations, particularly at higher floe counts.


# Acknowledgements

Our work is supported by the Office of Naval Research (ONR) grant
N00014-19-1-2421. The authors thank the authors of Subzero, Georgy Manucharyan
and Brandon Montemuro for their guidance and advice during the porting process.

# References

