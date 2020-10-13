import PyPlot.plot

"""
    plot_field(field::Array{<:Real,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap="coolwarm", vmin=nothing, vmax=nothing)

Plot the real-valued `field` on a unit cell with lattice vectors `a1` and `a2`.
"""
function plot_field(field::Array{<:Real,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap="coolwarm", vmin=nothing, vmax=nothing, label=nothing)
    ps = range(-0.5, stop=0.5, length=size(field,1)+1)
    qs = range(-0.5, stop=0.5, length=size(field,2)+1)
    xs = [p*a1[1]+q*a2[1] for p in ps, q in qs]
    ys = [p*a1[2]+q*a2[2] for p in ps, q in qs]
    if cmap == "coolwarm" && vmin == nothing && vmax == nothing
        vlim = maximum(abs.(field))
        vmin = -vlim
        vmax = +vlim
    end
    if cmap == "viridis" && vmin == nothing && vmax == nothing
        vmin = 0
        vmax = maximum(abs.(field))
    end
    pcolormesh(xs, ys, field, cmap=cmap, vmin=vmin, vmax=vmax)
    gca().set_aspect("equal")
    axis("off")
    colorbar()
    if label != nothing
        title(label)
    end
end

"""
    plot_field(field::Array{<:Real,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap="coolwarm", vmin=nothing, vmax=nothing)

Plot the complex-valued `field` on a unit cell with lattice vectors `a1` and `a2`.
"""
function plot_field(field::Array{<:Complex,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap="coolwarm", vmin=nothing, vmax=nothing, label=nothing)
    # Plot real and imaginary parts with the same color scale
    if cmap == "coolwarm" && vmin == nothing && vmax == nothing
        vlim = maximum(abs.(field))
        vmin = -vlim
        vmax = +vlim
    end
    if cmap == "viridis" && vmin == nothing && vmax == nothing
        vmin = 0
        vmax = maximum(abs.(field))
    end
    figure(figsize=(4,2))
    subplot(1,2,1)
    plot_field(real(field), a1, a2, cmap=cmap, vmin=vmin, vmax=vmax, label="Re[$label]")
    subplot(1,2,2)
    plot_field(imag(field), a1, a2, cmap=cmap, vmin=vmin, vmax=vmax, label="Im[$label]")
end


"""
    get_field(data::AbstractVector{<:Complex}, basis::PlaneWaveBasis;
                        k0=[0,0], t1s=-0.5:0.01:0.5, t2s=-0.5:0.01:0.5)

Convert the `data` from a `PlaneWaveBasis` to a real space grid.
"""
function get_field(data::AbstractVector{<:Complex}, basis::PlaneWaveBasis;
                    k0=[0,0], t1s=-0.5:0.01:0.5, t2s=-0.5:0.01:0.5)
    # Generate real space coordinates
    b1, b2 = basis.b1, basis.b2
    a1, a2 = bs_to_as(b1, b2)
    xs = [t1*a1[1]+t2*a2[1] for t1 in t1s, t2 in t2s]
    ys = [t1*a1[2]+t2*a2[2] for t1 in t1s, t2 in t2s]
    out = zeros(ComplexF64, length(t1s), length(t2s))
    # Inverse Fourier transform
    kxs = k0[1] .+ basis.kxs
    kys = k0[2] .+ basis.kys
    for (kx,ky,u) in zip(kxs, kys, data)
        for c in CartesianIndices(out)
            phase = kx*xs[c] + ky*ys[c]
            out[c] += u * exp(1im*phase)
        end
    end
    return out
end


"""
    plot_field(data::AbstractVector{<:Complex}, basis::PlaneWaveBasis;
                        k0=[0,0], cmap="coolwarm", vmin=nothing, vmax=nothing, label=nothing)

Converts the plane-wave amplitudes to real space using [`get_field`](@ref)
and then plots the field.
"""
function plot_field(data::AbstractVector{<:Complex}, basis::PlaneWaveBasis; k0=[0,0], cmap="coolwarm", vmin=nothing, vmax=nothing, label=nothing)
    field = get_field(data, basis, k0=k0)
    a1, a2 = bs_to_as(basis.b1, basis.b2)
    plot_field(field, a1, a2, cmap=cmap, vmin=vmin, vmax=vmax, label=label)
end


"""
    plot(geometry::Geometry)

Plot the permittivity and permeability of the `geometry` in real space.
"""
function PyPlot.plot(geometry::Geometry)
    # Plot permittivity
    plot_field(geometry.ep, geometry.a1, geometry.a2, cmap="viridis", label="ϵ")
    # Plot permeability
    plot_field(geometry.mu, geometry.a1, geometry.a2, cmap="viridis", label="μ")
end


"""
    plot(solver::Solver)

Plot the representation of the `solver`'s geometry in real space.

The `solver` approximates the geometry using a truncated basis of plane waves,
so `plot(solver)` lets us judge how accurately the geometry is represented.
"""
function PyPlot.plot(solver::Solver)
    col = findfirst((solver.basis.ps .==0) .& (solver.basis.qs .== 0))
    # Plot permittivity
    plot_field(solver.epc[:,col], solver.basis, cmap="viridis", label="ϵ")
    # Plot permeability
    plot_field(solver.muc[:,col], solver.basis, cmap="viridis", label="μ")
end


"""
    plot(mode::Eigenmode, [bloch_phase=true])

Plot the mode in real space.

To plot only the cell-periodic part of the Bloch wave, set `bloch_phase=false`.
"""
function PyPlot.plot(mode::Eigenmode; bloch_phase=true)
    if bloch_phase
        label = mode.data_label
        plot_field(mode.data, mode.basis, label=label, k0=mode.k0)
    else
        label = L"e^{-i k \cdot r} \cdot" * mode.data_label
        plot_field(mode.data, mode.basis, label=label)
    end
end
