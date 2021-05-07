![](assets/banner.png)

# Introduction

## Background

`Peacock.jl` - or the **P**lane-wave **E**xpansion **A**pproach to **C**haracterising **O**ptical **C**rystals in **k**-space - is a Julia package for studying photonic crystals using the Plane Wave Expansion Method.

Photonic crystals are materials whose optical properties arise from the structuring of the material when the size of the structures are comparable to the wavelengths of light. `Peacock.jl` is named for the irridescent colours of peacock feathers which arise not from pigmentation but from their photonic crystal structure, as shown below.

> ![Different zooms of a Peacock](assets/peacock_feathers_zoom.png)
> Peacock by allanlau2000 from pixabay. Feather by suju from pixabay. Electron microscope image of photonic crystal structure from Zi, Jian, et al. "Coloration strategies in peacock feathers.",  *Proceedings of the National Academy of Sciences* 100.22 (2003): 12576-12578. Copyright (2003) National Academy of Sciences.

As well as occuring naturally as in animals such as peacocks, advances in nanofabrication mean that 'designer' photonic crystals can now be manufactured for unprecedented control over the flow of light, with applications ranging from optical fibers to photonic circuitry. Photonic crystals are also a promising platform for more exotic materials like topological insulators.


## Features

Solve for...
* Transverse electric (TE) and transverse magnetic (TM) modes of 2D photonic crystals
  * Non-orthogonal unit cells
  * Inhomogeneous permittivity and/or permeability
* Chern numbers of topological photonic crystals using [built-in Wilson loop methods](@ref how_to_topology)


Focused on ease of use
* [Install](@ref getting_started_installation) with one line in Julia's package manager
* Simple visualisation of [geometry](@ref getting_started_geometry), [fields](@ref getting_started_modes), and [fully labelled band diagrams](@ref getting_started_bands)
* Reproduce and extend existing photonic crystals in the [`Peacock.Zoo`](@ref how_to_zoo) submodule
* Easily [accelerate calculations on CUDA-compatible GPUs](@ref how_to_GPU)


## Limitations

* Currently only implemented for 2D photonic crystals
* Like all methods that solve Maxwell's equations in Fourier space, the Plane Wave Expansion Method converges slowly for high contrast materials such as metals (``\epsilon < 0``)




## Contributors

* @sp94 (core)
* @kabume (GPU support)


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
