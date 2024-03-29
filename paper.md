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
    affiliation: 1
  - name: Mukund Gupta
    orcid: 0000-0000-0000-0000
    affiliation: 2
  - name: Andrew Thomspn
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: California Institute of Technology, CA, USA
   index: 1
 - name: Delft Institute of Technology
   index: 2
date: 1 April 2024
bibliography: paper.bib
---

# Summary

Arctic sea ice extent and concentration continue to decline at rates that are commonly underestimated by climate projection models. Potential sources of uncertainty arise from an inaccurate representation of interactions between the ocean and sea ice within these climate models, as well as the simplification of sea-ice dynamics for the sake of reducing computational complexity. Discrete-element models (DEMs), where each piece of sea ice is represented as an individual simulation element, all of which can dynamically interact, are used to explore these uncertainties and study fine-scale sea ice dynamics. However, these models are computationally expensive and have not yet been coupled to a dynamic ocean to explore two-way feedbacks.

# Statement of need

SubZero [@Montemuro:2023], a novel DEM written in MATLAB, pushes beyond
traditional models by representing sea ice floes as polygonal elements that
change in shape, mass, and number over time as a result of interactions with
other floes and topographical elements. These features address the uncertainty
of simplified sea ice dynamics within continuous models mentioned above, but are
computationally challenging. To increase the scale and speed of simulations, it
was determined that SubZero should be ported from MATLAB to the Julia
programming language [@Bezanson:2017] and re-engineered to improve both
performance and usability. Additionally, exploration of sea ice and ocean
dynamics, the second source of uncertiantly discussed above, requires coupling
with a dynamic ocean model. Porting SubZero to Julia and adding two-way coupling
infrastructure allows coupling with Oceananigans.jl [@Ramadhan:2020], creating a
way to explore new scientific questions.

Subzero.jl, the redesigned version of SubZero in Julia, achieves speeds up to
[NEED TO KNOW] times faster than the original model when running the Nare's
Straight simulation discussed in detail in @Manucharyan:2022 and
@Montemuro:2023. This speedup will enable longer runs, allowing study of an
entire year of sea ice evolution within [] hours.
[ something something performance benchmarks ] 

Furthermore, the new coupling framework allows...

Finally, the switch to Julia allows a more extendable, user-friendly interface.
With a modular simulation and model object, users can craft detailed simulations
with a script-based interface without ever interacting with the source code.
Furthermore, it is easy for the user to extend simulation functions, such as the
creation of new domain boundary types and floe fracture criteria, within
their own scripts, rather than within the source code, using Julia's multiple
dispatch paradigm.

Subzero.jl is continuously tested against a suite of unit tests and integration tests that
compare its behavior to the original MATLAB model, and confirm that the model conserves both
energy and momentum. 

# Citations


# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

Our work is supported by the Office of Naval Research (ONR) grant
N00014-19-1-2421. The authors thank the authors of Subzero, Georgy Manucharyan
and Brandon Montemuro for their guidance and advice during the porting process.

# References