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


# Example usage

![Example reproduction of a fragilely topological photonic crystal
[@blanco2019engineering]. (a) Unit cell of the photonic crystal.
(b-d) Out-of-plane component of the electric field of the first three
transverse-magnetic polarised modes at $\Gamma$. (e) Band diagram of the
transverse-magnetic polarised modes. (f) The Wilson loop spectrum of bands 2-3
wind, indicating non-trivial band topology. (g) The Wilson loop spectrum of the
full valence band space does not wind, indicating that bands 2-3 are 'fragilely'
topological. \label{fig:examples}](figures/examples.pdf)

## Loading from the `Peacock.Zoo` submodule

In this example we reproduce results that demonstrate the photonic crystal with
'fragile' band topology that was introduced in [@blanco2019engineering].

This crystal can be loaded from the `Peacock.Zoo` submodule using `make_dePaz_frag`.
```julia
using Peacock, Peacock.Zoo, Parameters

# Size of the plane-wave basis
fourier_space_cutoff = 7

# Load the fragile photonic topological insulator (Blanco de Paz et al 2019)
@unpack geometry, solver, polarisation = 
        make_dePaz_frag(fourier_space_cutoff)

# Visualise the geometry
plot(geometry)
```

If you model your own photonic crystal with `Peacock.jl`, you can add your geometry
to the `Zoo` submodule to help others reproduce your work.


## Plotting the Wilson loop winding

A winding in the Wilson loop spectrum can indicate a non-trivial topological
phase [@blanco2020tutorial], with the Chern number given by the winding. 

First, we define the ``k``-path we want to scan along, labelling the high
symmetry points using `BrillouinZoneCoordinate`.
```julia
# The Wilson loops are (by default) along b2, so we define a straight
# path from Γ to Γ+b1 - we will scan along this path
ks = [
    BrillouinZoneCoordinate(0.0, 0.0, "Γ"),
    BrillouinZoneCoordinate(0.5, 0.0, "M"),
    BrillouinZoneCoordinate(1.0, 0.0, "Γ")
]
```

Now we can reproduce the Wilson loop winding of \autoref{fig:examples}f-g.
```julia
# Wilson loop of all three valence bands
figure(figsize=(3,2))
plot_wilson_loop_winding(solver, ks, polarisation, 1:3, dk_outer=0.25)
title("Bands 1-3")

# Wilson loop of just the second and third bands
figure(figsize=(3,2))
plot_wilson_loop_winding(solver, ks, polarisation, 2:3, dk_outer=0.25)
title("Bands 2&3")
```

In the first figure, the Wilson loops through the Hilbert spaces of bands 2&3
wind with Chern numbers ±1, indicating some non-trivial topology.
However, the second figure shows that including the (trivial) acoustic band in
the Wilson loop calculation removes the topological winding, and consequently
bands 2&3 are said to be 'fragilely topological' [@blanco2019engineering].


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