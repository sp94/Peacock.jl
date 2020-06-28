<div align="center">
<img src="docs/src/assets/banner.png" alt="Peacock logo"></img>
</div>

---

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sp94.github.io/Peacock.jl/dev)

## Introduction

`Peacock.jl` - or the **P**lane-wave **E**xpansion **A**pproach to **C**haracterising **O**ptical **C**rystals in **k**-space - is a Julia package for studying photonic crystals using the Plane Wave Expansion Method.

Photonic crystals are materials whose optical properties arise from the structuring of the material when the size of the structures are comparable to the wavelengths of light. `Peacock.jl` is named for the irridescent colours of peacock feathers which arise not from pigmentation but from their photonic crystal structure, as shown below.

> ![Different zooms of a Peacock](docs/src/assets/peacock_feathers_zoom.png)
> Peacock by allanlau2000 from pixabay. Feather by suju from pixabay. Electron microscope image of photonic crystal structure from Zi, Jian, et al. "Coloration strategies in peacock feathers.",  *Proceedings of the National Academy of Sciences* 100.22 (2003): 12576-12578. Copyright (2003) National Academy of Sciences.

As well as occuring naturally as in animals such as peacocks, advances in nanofabrication mean that 'designer' photonic crystals can now be manufactured for unprecedented control over the flow of light, with applications ranging from optical fibers to photonic circuitry. Photonic crystals are also a promising platform for more exotic materials like topological insulators.


## Features

Solve for...
* Transverse electric (TE) and transverse magnetic (TM) modes of 2D photonic crystals
  * Non-orthogonal unit cells
  * Inhomogeneous permittivity and/or permeability
* Chern numbers of topological photonic crystals using [built-in Wilson loop methods](https://sp94.github.io/Peacock.jl/dev/how-tos/wilson_loops)


Focused on ease of use
* [Install](@ref https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_installation-1) with one line in Julia's package manager
* Simple visualisation of [geometry](@ref https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_geometry-1), [fields](@ref https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_modes-1), and [fully labelled band diagrams](@ref https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_bands-1)
* Reproduce and extend existing photonic crystals in the [`Peacock.Zoo`](@ref https://sp94.github.io/Peacock.jl/dev/how-tos/zoo/#how_to_zoo-1) submodule


## Limitations

* Currently only implemented for 2D photonic crystals
* Like all methods that solve Maxwell's equations in Fourier space, the Plane Wave Expansion Method converges slowly for high contrast materials such as metals (``\epsilon < 0``)
