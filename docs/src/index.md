---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "Subzero.jl"
  text: ""
  tagline: Fast and Flexible Sea Ice Dynamics
  actions:
    - theme: brand
      text: Introduction
      link: /introduction
    - theme: alt
      text: View on Github
      link: https://github.com/Caltech-OCTO/Subzero.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - icon: <img width="64" height="64" src="https://rawcdn.githack.com/JuliaLang/julia-logo-graphics/f3a09eb033b653970c5b8412e7755e3c7d78db9e/images/juliadots.iconset/icon_512x512.png" alt="Julia code"/>
    title: Pure Julia code
    details: Fast, understandable, extensible functions
    link: /introduction
  - icon: <img width="64" height="64" src="https://fredrikekre.github.io/Literate.jl/v2/assets/logo.png" />
    title: Literate programming
    details: Documented source code with examples!
    link: /source/methods/clipping/cut
---


<p style="margin-bottom:2cm"></p>

<div class="vp-doc" style="width:80%; margin:auto">

# What is Subzero.jl?

Subzero is an easy-to-use Julia translation of Manucharyan and Montemuro’s model described in the paper “SubZero: A Sea Ice Model With an Explicit Representation of the Floe Life Cycle.” The model has been restructured to leverage Julia’s language abstractions for ease of setting up new simulation runs and allowing more types of simulations without code modifications. It has been designed for stand-alone simulations.

Subzero.jl was ported and restructured by Skylar Gering and originally developed by Georgy Manucharyan and Brandon Montemuro.

We welcome contributions, either as pull requests or discussion on issues!

</div>
