---
title: '`Peacock.jl`: Photonic crystals in Julia'
tags:
  - Julia
  - physics
  - optics
  - photonics
  - photonic crystals
  - photonic topological insulators
  - Maxwell's equations
  - metamaterials
  - plane wave expansion
authors:
  - name: Samuel J. Palmer
    orcid: 0000-0003-0485-7047
    affiliation: "1"
  - name: Vincenzo Giannini
    orcid: 0000-0001-8025-4964
    affiliation: "2"
affiliations:
 - name: The Blackett Laboratory, Imperial College London, London, SW7 2AZ, UK
   index: 1
 - name: Instituto de Estructura de la Materia (IEM), Consejo Superior de Investigaciones Científicas (CSIC), Serrano 121, 28006, Madrid, Spain
   index: 2
date: 21 May 2020
bibliography: paper.bib
---

# Summary

![The irridescent colours of peacocks arise from the nanoscale 'photonic
crystal' structure of the feathers, rather than from pigmentation. The image of
the peacock's photonic crystal structure is reproduced from @zi2003coloration,
Copyright National Academy of Sciences. \label{fig:zoom}](../docs/src/assets/peacock_feathers_zoom.png)

The **P**lane-wave **E**xpansion **A**pproach to **C**haracterising **O**ptical
**C**rystals in **k**-space - otherwise known as
`Peacock.jl` - is a Julia package for studying
photonic crystals using the Plane Wave Expansion Method [@rumpf2006design].
A photonic crystal is a material whose optical properties arise from its
periodic structure [@yablonovitch1987inhibited; @john1987strong], and
`Peacock.jl` is named for the irridescent colours of peacock feathers which
arise not from pigmentation but from the scattering of light by photonic crystals
[@zi2003coloration], as shown in \autoref{fig:zoom}.

The response of a photonic crystal is strongest
when the periodicity of the structure is comparable to the wavelength of light.
For visible light, photonic crystals are built from components that are just a
few hundred nanometers in size. Advances in nanofabrication mean that 'designer'
photonic crystals can now be manufactured for unprecedented control over the
flow of light, with applications ranging from optical fibers
[@knight2003photonic] to photonic circuitry [@joannopoulos2008molding].
Photonic crystals are also a promising platform for materials known as
photonic topological insulators [@wu2015scheme; @wang2019band;
@blanco2019engineering]. These are the photonic analogue of the electronic topological
insulators [@kane2005z; @kane2005quantum; @bernevig2006quantum]
for which the 2016 Nobel Prize in Physics was awarded, and may allow
light to be guided around defects, impurities, and sharp corners without
backscattering [@rider2019perspective].

# Statement of need

`Peacock.jl` provides a user-friendly interface to calculate and analyse
the eigensolutions of 2D photonic crystals,
with support for non-orthogonal unit cells and inhomogeneous permittivity and/or
permeability. As well as the common tools for eigenmode analysis,
such as visualising the eigenmodes or the plotting the eigenvalues on band diagrams as
in \autoref{fig:examples}a-e, `Peacock.jl` also includes built-in using Wilson loop
methods [@blanco2020tutorial] to study the topology of photonic bands,
as in \autoref{fig:examples}f-g.
Although there exists open-source software to study photonic bands [@johnson2001block]
and to study topology in electronic bands [@gresch2017z2pack], to our knowledge
`Peacock.jl` is the first open-source package for studying the
band topology in photonic crystals.

# Examples

Finally, `Peacock.jl` includes a `Zoo` submodule of crystals from published
works. If you model a photonic crystal with `Peacock.jl`, please consider
contributing your geometry to the `Zoo` submodule to help others reproduce your work.

![Example reproduction of a fragilely topological photonic crystal
[@blanco2019engineering]. (a) Unit cell of the photonic crystal.
(b-d) Out-of-plane component of the electric field of the first three
transverse-magnetic polarised modes at $\Gamma$. (e) Band diagram of the
transverse-magnetic polarised modes. (f) The Wilson loop spectrum of bands 2-3
wind, indicating non-trivial band topology. (g) The Wilson loop spectrum of the
full valence band space does not wind, indicating that bands 2-3 are 'fragilely'
topological. \label{fig:examples}](examples.pdf)

# Acknowledgements

S.J.P. acknowledges his studentship from the Centre for Doctoral Training on
Theory and Simulation of Materials at Imperial College London funded
by EPSRC Grant No. EP/L015579/1.
​
V.G. acknowledges the Spanish Ministerio de Economia y Competitividad for
financial support through the grant NANOTOPO (FIS2017-91413-EXP) and also the
Ministerio de Ciencia, Innovació n y Universidades through the grant MELODIA
(PGC2018-095777-B-C21).

# References