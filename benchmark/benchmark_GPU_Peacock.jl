using PyPlot
using Peacock
using CUDA
GPU = 1

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

figure(figsize=(4,3))
@time plot_band_diagram(solver, ks, TE, GPU, color="red",
            bands=1:4, dk=0.1, frequency_scale=1/2pi)
@time plot_band_diagram(solver, ks, TM, GPU, color="blue",
            bands=1:4, dk=0.1, frequency_scale=1/2pi)
ylim(0,0.8)
