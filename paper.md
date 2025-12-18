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
    orcid: 0000-0000-0000-0000
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

# Statement of need
Arctic sea ice extent and concentration continue to decline at rates that are
commonly underestimated by climate projection models. Potential sources of
uncertainty arise from an inaccurate representation of interactions between the
ocean and sea ice within these climate models, as well as the simplification of
sea-ice dynamics for the sake of reducing computational complexity.
Discrete-element models (DEMs), where each piece of sea ice is represented as an
individual simulation element, all of which can dynamically interact, are used
to explore these uncertainties and study fine-scale sea ice dynamics. However,
these models are computationally expensive and are not commonly coupled to a
dynamic ocean to explore two-way feedbacks. We present Subzero.jl, a native Julia [@Bezanson:2017]
version of MATLAB discrete-element model SubZero [@Montemuro:2023] that addresses both of these two
problems. 

# Summary
SubZero [@Montemuro:2023], a novel DEM written in MATLAB, pushes beyond
traditional models by representing sea ice floes as polygonal elements that
change in shape, mass, and number over time as a result of interactions with
other floes and topographical elements. These features address the uncertainty
of simplified sea ice dynamics within continuous models mentioned above, but are
computationally expensive. To increase the scale and speed of simulations, it
was determined that SubZero should be ported from MATLAB to the Julia
programming language [@Bezanson:2017] and re-engineered to improve both
performance and usability. Additionally, exploration of sea ice and ocean
dynamics, the second source of uncertiantly discussed above, requires coupling
with a dynamic ocean model. Porting SubZero to Julia and adding two-way coupling
infrastructure allows coupling with Oceananigans.jl [@Ramadhan:2020], allowing
exploration of new scientific questions.

First of all, the switch to Julia allows a more extendable, user-friendly interface.
With a modular simulation and model object, Subzero.jl allows users to craft detailed
simulations with a script-based interface without ever interacting with the source code.
Users can easily turn on and off various physical processes like fracturing, ridging, and welding.
It is also easy for the user to extend existing functionality, such as the creating new domain boundary types and floe fracture criteria, within their own scripts, rather than within the source code, using Julia's multiple dispatch paradigm. This flexibility is highlighted in the tutorial and examples
provided with the documentation.

Subzero.jl also achieves speeds up to [ADD DETAILS] times faster than the original model when running the repository example scripts. This speed up will enable longer runs, allowing ...

Furthermore, the new coupling framework allows...

Finally, Subzero.jl is continuously tested against a suite of unit tests and integration
tests that compare its behavior to the original MATLAB model, and confirm that
the model conserves both energy and momentum. 


<!-- Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

Our work is supported by the Office of Naval Research (ONR) grant
N00014-19-1-2421. The authors thank the authors of Subzero, Georgy Manucharyan
and Brandon Montemuro for their guidance and advice during the porting process.

# References