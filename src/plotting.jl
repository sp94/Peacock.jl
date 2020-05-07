import PyPlot.plot

"""
To do
"""
function plot_field(field::Array{<:Real,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap="coolwarm", vmin=nothing, vmax=nothing)
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
end


"""
To do
"""
function plot_field(field::Array{<:Complex,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap="coolwarm", vmin=nothing, vmax=nothing)
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
    plot_field(real(field), a1, a2, cmap=cmap, vmin=vmin, vmax=vmax)
    subplot(1,2,2)
    plot_field(imag(field), a1, a2, cmap=cmap, vmin=vmin, vmax=vmax)
end


"""
To do
"""
function plot_field(u::AbstractVector{<:Complex}, basis::PlaneWaveBasis; cmap="coolwarm", vmin=nothing, vmax=nothing, k0=[0,0])
    field = get_field(u, basis, k0=k0)
    a1, a2 = bs_to_as(basis.b1, basis.b2)
    plot_field(field, a1, a2, cmap=cmap, vmin=vmin, vmax=vmax)
end


"""
To do
"""
function PyPlot.plot(geometry::Geometry)
    # Plot permittivity
    plot_field(geometry.ep, geometry.a1, geometry.a2, cmap="viridis")
    subplot(1,2,1)
    title("Re[ϵ]")
    colorbar()
    subplot(1,2,2)
    title("Im[ϵ]")
    colorbar()
    # Plot permeability
    plot_field(geometry.mu, geometry.a1, geometry.a2, cmap="viridis")
    subplot(1,2,1)
    title("Re[μ]")
    colorbar()
    subplot(1,2,2)
    title("Im[μ]")
    colorbar()
end


"""
To do
"""
function PyPlot.plot(solver::Solver)
    col = findfirst((solver.basis.ps .==0) .& (solver.basis.qs .== 0))
    # Plot permittivity
    plot_field(solver.epc[:,col], solver.basis, cmap="viridis")
    subplot(1,2,1)
    title("Re[ϵ]")
    colorbar()
    subplot(1,2,2)
    title("Im[ϵ]")
    colorbar()
    figure()
    # Plot permeability
    plot_field(solver.muc[:,col], solver.basis, cmap="viridis")
    subplot(1,2,1)
    title("Re[μ]")
    colorbar()
    subplot(1,2,2)
    title("Im[μ]")
    colorbar()
end


"""
To do
"""
function PyPlot.plot(mode::Mode; bloch_phase=true)
    if bloch_phase
        plot_field(mode.data, mode.basis, k0=mode.k0)
    else
        plot_field(mode.data, mode.basis)
    end
    if bloch_phase
        label = mode.label
    else
        label = L"e^{-i k \cdot r} \cdot" * mode.label
    end
    subplot(1,2,1)
    title("Re[$label]")
    colorbar()
    subplot(1,2,2)
    title("Im[$label]")
    colorbar()
end
