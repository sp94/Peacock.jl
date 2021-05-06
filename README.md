<div align="center">
<img src="docs/src/assets/banner.png" alt="Peacock logo"></img>
</div>

---

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sp94.github.io/Peacock.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sp94.github.io/Peacock.jl/dev)
[![Coverage Status](https://coveralls.io/repos/github/sp94/Peacock.jl/badge.svg?branch=master)](https://coveralls.io/github/sp94/Peacock.jl?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02678/status.svg)](https://doi.org/10.21105/joss.02678)

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
* [Install](https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_installation-1) with one line in Julia's package manager
* Simple visualisation of [geometry](https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_geometry-1), [fields](https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_modes-1), and [fully labelled band diagrams](https://sp94.github.io/Peacock.jl/dev/tutorials/getting_started/#getting_started_bands-1)
* Reproduce and extend existing photonic crystals in the [`Peacock.Zoo`](https://sp94.github.io/Peacock.jl/dev/how-tos/zoo/#how_to_zoo-1) submodule
* Set `use_GPU=true` for accelerated calculations on CUDA-compatible GPUs


## Limitations

* Currently only implemented for 2D photonic crystals
* Like all methods that solve Maxwell's equations in Fourier space, the Plane Wave Expansion Method converges slowly for high contrast materials such as metals (ϵ < 0)


## Referencing

If you use `Peacock.jl` in your work, please consider citing us as
```
@article{palmer2020peacock,
  title={Peacock.jl: Photonic crystals in {Julia}},
  author={Palmer, Samuel J and Giannini, Vincenzo},
  journal={Journal of Open Source Software},
  volume={5},
  number={54},
  pages={2678},
  year={2020}
}
```
