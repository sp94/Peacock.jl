using PyPlot

function plot_field(field::Array{<:Real,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap=nothing, vmin=nothing, vmax=nothing)
    ps = range(-0.5, stop=0.5, length=size(field,1))
    qs = range(-0.5, stop=0.5, length=size(field,2))
    xs = [p*a1[1]+q*a2[1] for p in ps, q in qs]
    ys = [p*a1[2]+q*a2[2] for p in ps, q in qs]
    if cmap == "coolwarm" && vmin == nothing && vmax == nothing
        vlim = maximum(abs.(field))
        vmin = -vlim
        vmax = +vlim
    end
    pcolormesh(xs, ys, field, cmap=cmap, vmin=vmin, vmax=vmax)
    gca().set_aspect("equal")
end

function plot_field(field::Array{<:Complex,2}, a1::Array{<:Real,1}, a2::Array{<:Real,1}; cmap=nothing, vmin=nothing, vmax=nothing)
    # Plot real and imaginary parts with the same color scale
    if cmap == "coolwarm" && vmin == nothing && vmax == nothing
        vlim = maximum(abs.(field))
        vmin = -vlim
        vmax = +vlim
    end
    subplot(1,2,1)
    plot_field(field, a1, a2, cmap=cmap, vmin=vmin, vmax=vmax)
    subplot(1,2,2)
    plot_field(field, a1, a2, cmap=cmap, vmin=vmin, vmax=vmax)
end

"""
Plot a plane wave expansion as a real space field.
"""
function plot_mode(cg::ConvolvedGeometry, u::Array{<:Complex,1}, k::Array{<:Real,1})
    field = get_field(cg, u, k0=k)
    plot_field(field, cg.a1, cg.a2, cmap="coolwarm")
end
