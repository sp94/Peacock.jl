### Issue
The computer cost will very large if we extent Peacock.jl to 3D or need a large cutoff. Using parallel technology can significantly reduce the computing time. Can we parallelize the solver?

### GPU parallel
First,  We can use [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) to make computing run on Nvidia GPU. [CuSolver](https://docs.nvidia.com/cuda/cusolver/index.html) is a linear algebra library of CUDA. It is easy to solve dense linear systems of the form `Ax = b`. We only need to convert Julia Array to CuArray, almost all function of LinearAlgebra.jl can be running on GPU:
```julia
# running on CPU
a = rand(2, 2); b = rand(2, 2)
c = a / b
# running on GPU
using CUDA
a_d = CuArray(a); b_d = CuArray(b)
c_d = a_d / b_d
```
Unfortunately, it may be inefficient to use CUDA to solve the eigenvalue and eigenvector of dense matrix ([discourse1](https://discourse.julialang.org/t/computing-eigenvalues-eigenvectors-using-gpu/12396) [discourse2](https://discourse.julialang.org/t/eigenvalues-for-lots-of-small-matrices-gpu-batched-vs-cpu-eigen/50792)). I also run  benchmarks :
```
julia>  a= rand(100,100); b = rand(100,100);
julia> @btime CUDA.CUSOLVER.syevjBatched!('V','U',CuArray(b)/CuArray(a));
  36.093 ms (577 allocations: 11.58 KiB)

julia> @btime eigen(a,b);
  29.361 ms (21 allocations: 428.78 KiB)
```

It should be note that the function `CUDA.eigen` of CUDA.jl only can run with Julia Array and just same as `LinearAlgebra.eigen` of LinearAlgebra.jl. By the wave, `syevjBatched!` only support `Float32/Float64`. 

But we can still speed up Peacock.jl by using CUDA.jl, since GPU can solve `Ax = b` efficiently. I focus on the model of [square lattice](https://sp94.github.io/Peacock.jl/stable/tutorials/getting_started/): 

```julia
using PyPlot
using Peacock
using LinearAlgebra
using SparseArrays
using BenchmarkTools
using CUDA

get_k = Peacock.get_k
DiagonalMatrix = Peacock.DiagonalMatrix

# Permittivity
function epf(x,y)
    # equation of a circle with radius 0.2a
    if x^2+y^2 <= 0.2^2
        # dielectric inside the circle
        return 8.9
    else
        # air outside the circle
        return 1
    end
end

# Permeability is unity everywhere
function muf(x,y)
    return 1
end

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector
d1 = 0.01  # resolution along first lattice vector
d2 = 0.01  # resolution along second lattice vector
geometry = Geometry(epf, muf, a1, a2, d1, d2)
fourier_space_cutoff = 17 #greater than tutorials
solver = Solver(geometry, fourier_space_cutoff)
G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]
ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]

polarisation = TE
k = ks[1]
basis = solver.basis
epc, muc = solver.epc, solver.muc
Kx = DiagonalMatrix(basis.kxs) + k[1]*I
Ky = DiagonalMatrix(basis.kys) + k[2]*I
Kx_d = CuArray(Kx); Ky_d = CuArray(Ky);
epc_d = CuArray(epc); muc_d = CuArray(muc)

```

and run benchmarks like that:

```
julia> @btime LHS_d = Kx_d/epc_d*Kx_d + Ky_d/epc_d*Ky_d;
  22.768 ms (20762 allocations: 337.38 KiB)

julia> @btime LHS_d = Kx_d/epc_d*Kx_d + Ky_d/epc_d*Ky_d;
  22.931 ms (18741 allocations: 303.83 KiB)

julia> @btime LHS = Kx/epc*Kx + Ky/epc*Ky;
  2.433 s (20 allocations: 6.96 MiB)

julia> @btime LHS = Kx/epc*Kx + Ky/epc*Ky;
  2.660 s (20 allocations: 6.96 MiB)

julia> size(LHS)
(225, 225)
```

GPU-based program is about 90 times faster than CPU! The other benchmark is test `eigen`:

```
julia> @btime eigen(RHS \ LHS);
  1.362 s (29 allocations: 3.55 MiB)

julia> @btime eigen(LHS, RHS);
  542.625 ms (19 allocations: 2.46 MiB)

julia> @time eigen(Array(RHS_d \ LHS_d));
  0.209265 seconds (17.91 k allocations: 3.053 MiB)
```

Despite `eigen(LHS, RHS)` is more efficient  than `eigen(RHS \ LHS)`,  `eigen(Array(RHS_d \ LHS_d))` is still the fastest one! 

Now, we can change the code of `solver.jl`. Only one function is added:

```julia
function solve(solver::Solver, k::AbstractVector{<:Real}, polarisation::Polarisation, GPU::Int; bands=:)
    basis = solver.basis
    if GPU == 1
        epc, muc = CuArray(solver.epc), CuArray(solver.muc)
        Kx = CuArray(DiagonalMatrix(basis.kxs) + k[1]*I)
        Ky = CuArray(DiagonalMatrix(basis.kys) + k[2]*I)
        if polarisation == TE
            LHS = Kx/epc*Kx + Ky/epc*Ky
            RHS = muc
            label = L"H_z"
        elseif polarisation == TM
            LHS = Kx/muc*Kx + Ky/muc*Ky
            RHS = epc
            label = L"E_z"
        end
        freqs_squared, modes_data = eigen(Array(RHS \ LHS));
    else
        epc, muc = solver.epc, solver.muc
        Kx = DiagonalMatrix(basis.kxs) + k[1]*I
        Ky = DiagonalMatrix(basis.kys) + k[2]*I
        if polarisation == TE
            LHS = Kx/epc*Kx + Ky/epc*Ky
            RHS = muc
            label = L"H_z"
        elseif polarisation == TM
            LHS = Kx/muc*Kx + Ky/muc*Ky
            RHS = epc
            label = L"E_z"
        end
        freqs_squared, modes_data = try
            eigen(LHS, RHS)
        catch
            eigen(RHS \ LHS)
        end
    end
    
    freqs = sqrt.(freqs_squared)
    # Sort by increasing frequency
    idx = sortperm(freqs, by=real)
    freqs = freqs[idx][bands]
    modes_data = modes_data[:,idx][:,bands]
    # Eigenmodes are weighted by the RHS of the generalised eigenvalue problem
    weighting = RHS
    modes = Mode[]
    for i in 1:length(freqs)
        mode = Mode(k, freqs[i], modes_data[:,i], weighting, basis, label)
        push!(modes, mode)
    end
    return modes
end
```

The last function of `band_diagrams.jl` is change to:

```julia
function plot_band_diagram(solver::Solver, ks, polarisation::Polarisation, GPU::Int;
            dk=nothing, labels=[], bands=1:10, frequency_scale=1, color="k", markersize=nothing)
    # Convert BrillouinZoneCoordinate to labelled positions in k space
    if labels == []
        labels = [hasproperty(x,:label) ? x.label : "" for x in ks]
    end
    ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]
    # Wrap all the variables into a single function of k that returns frequencies
    function my_solve(k)
        modes = solve(solver, k, polarisation, GPU)
        return [mode.frequency for mode in modes]
    end
    # Pass on to more general function
    plot_band_diagram(my_solve, ks; dk=dk, labels=labels, bands=bands, frequency_scale=frequency_scale, color=color, markersize=markersize)
end
```

And we can run the last benchmark -- plot a band diagram, the model is just same as before:

```julia
using PyPlot
using Peacock
using CUDA
GPU = 1 # running on GPU

function epf(x,y)
    # equation of a circle with radius 0.2a
    if x^2+y^2 <= 0.2^2
        # dielectric inside the circle
        return 8.9
    else
        # air outside the circle
        return 1
    end
end

# Permeability is unity everywhere
function muf(x,y)
    return 1
end

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector
d1 = 0.01  # resolution along first lattice vector
d2 = 0.01  # resolution along second lattice vector
geometry = Geometry(epf, muf, a1, a2, d1, d2)
fourier_space_cutoff = 7
solver = Solver(geometry, fourier_space_cutoff)
G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]
```

The results when running on terminal:

```
julia> @time plot_band_diagram(solver, ks, TE, 1, color="red",
                   bands=1:4, dk=0.1, frequency_scale=1/2pi)
  3.792055 seconds (510.63 k allocations: 124.970 MiB, 0.71% gc time, 0.14% compilation time)
(0, 0.9109956614327619)

julia> @time plot_band_diagram(solver, ks, TM, 1, color="red",
                   bands=1:4, dk=0.1, frequency_scale=1/2pi)
  3.321212 seconds (571.62 k allocations: 125.828 MiB, 0.57% gc time)
(0, 0.9109956614327619)

julia> @time plot_band_diagram(solver, ks, TE, 0, color="red",
                   bands=1:4, dk=0.1, frequency_scale=1/2pi)
 29.343061 seconds (296.76 k allocations: 50.585 MiB, 0.03% gc time)
(0, 0.9109956614327619)

julia> @time plot_band_diagram(solver, ks, TM, 0, color="red",
                   bands=1:4, dk=0.1, frequency_scale=1/2pi)
 28.964005 seconds (294.21 k allocations: 50.534 MiB, 0.03% gc time)
(0, 0.9109956614327619)
```

When running on terminal, the diagram will be displayed in real time at one point by one, which may reduce efficiency. So I also run this benchmarks on Jupyter:



![](https://github.com/kabume/Peacock.jl/blob/master/benchmark/GPU_0.png)

![](https://github.com/kabume/Peacock.jl/blob/master/benchmark/GPU_1.png)

GPU-based computing is about 30 times faster than CPU-based.

### CPU Parallel(TODO

We can also use Multi-Threading (`@threads`) or Multi-Programming (`using Distributed; @distributed`) to parallel the `for` loop [like](https://github.com/sp94/Peacock.jl/blob/master/src/band_diagrams.jl#L56):

```julia
for (x,k) in zip(xs,ks)
    frequencies = my_solve(k)
    xs_k = [x for _ in 1:length(frequencies[bands])]
    ys = frequency_scale * frequencies[bands]
    PyPlot.plot(xs_k, real(ys), ".", color=color, alpha=0.5, markersize=markersize)
end
```

