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
Subzero.jl is a discrete-element model (DEM) for simulating sea-ice floe dynamics, designed to be fast, flexible, and user-friendly.
The model is a port of the MATLAB-based SubZero [@Montemuro2023] to the Julia programming language [@Bezanson:2017], with substantial re-engineering to increase speed, to improve extensibility and usability, and to provide native infrastructure for explicit coupling with external ocean dynamics models.
The model is available as a registered Julia package.

# Statement of need
Arctic sea ice extent and concentration continue to decline at rates that are commonly underestimated by climate projection models [@cite].
Potential sources of uncertainty arise from an inaccurate representation of interactions between the ocean and sea ice within these climate models, as well as the simplification of sea-ice dynamics for the sake of reducing computational complexity.
Discrete-element models (DEMs), where each piece of sea ice is represented as an individual simulation element, all of which can dynamically interact, are used to explore these uncertainties and study fine-scale sea ice dynamics. 
While a range of sea ice DEMs are available [e.g., @Hopkins2004,@Herman2013,@Rabatel2015,@Damsgaard2018, and others], many of these models make a range of simplifications, particularly in the geometric representation of the floes.
The MATLAB-based SubZero model [@Montemuro2023] allows for complex, and evolving, floe shapes but is computationally expensive. 
Moreover, due to the highly connected nature of sea ice floes and ocean processes [@Horvat2018,@Gupta2022,@Gupta2024,@Brenner2023c], there is a need for DEM simulations to be coupled to a dynamic ocean model to explore two-way feedbacks--a feature not readily available for most extant models. 
We present Subzero.jl, a native Julia [@Bezanson2017] version of the MATLAB discrete-element model SubZero [@Montemuro2023, @Manucharyan2022] that addresses both of these two problems. 

# Functionality
Subzero.jl, represents sea ice floes as polygonal elements that can change in shape, mass, and number over time as a result of interactions with other floes and topographical elements, and move in response to forcing from either the atmosphere or the ocean [@Manucharyan2022].

This new Julia implementation allows a more extendable, user-friendly interface relative to the previous MATLAB version.
With a modular simulation and model object, Subzero.jl allows users to craft detailed simulations with a script-based interface without ever interacting with the source code.
Users can easily turn on and off various physical processes, such as fracturing, ridging, and welding.
It is also easy for the user to extend existing functionality, such as creating new domain boundary types and floe fracture criteria, within their own scripts, rather than within the source code, using Julia's multiple dispatch paradigm. 
This flexibility is highlighted in the tutorial and examples provided with the documentation.

Subzero.jl also achieves speeds up to [ADD DETAILS] times faster than the original model when running the repository example scripts (see \autoref{fig:speed_comparison}). 
This speed up will enable longer and more complex simulations, allowing ...

![Caption for example figure.\label{fig:speed_comparison}](speed_comparison.png)

Furthermore, the new coupling framework enables two-way coupled simulations with the large-eddy simulation (LES) ocean model Oceananigans.jl [@Ramadhan2020,@Wagner2025], enabling exploration of new scientific questions.
In this configuration, sea-ice floes in Subzero.jl are forced by gridded ocean velocity fields supplied by Oceananigans.jl, while floe-resolved ice-ocean stresses computed by Subzero.jl are re-gridded and returned to the ocean model as a spatially heterogeneous surface boundary condition.
The exchange of fields is implemented using callback functionality in Oceananigans.jl, allowing serial coupling at user-defined frequencies.
This coupled DEM-LES system resolves the mechanical interaction between evolving ocean fields and individual sea ice floes, enabling investigation of floe-scale modulation of upper-ocean dynamics [e.g., @Brenner2025].
An example coupled simulation is demonstrated in \autoref{fig:oceananigans_coupled_example}.

![Example Subzero-Oceananigans coupled simulation.
\label{fig:oceananigans_coupled_example}](oceananigans_coupled_example.png)

Finally, Subzero.jl is continuously tested against a suite of unit tests and integration tests that compare its behavior to the original MATLAB model, and confirm that the model conserves both energy and momentum.