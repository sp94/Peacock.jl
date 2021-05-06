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
fourier_space_cutoff = 17
solver = Solver(geometry, fourier_space_cutoff)
G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]
ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]

polarisation = TE

k = ks[1]
bands = 1:2

basis = solver.basis
epc, muc = solver.epc, solver.muc
Kx = DiagonalMatrix(basis.kxs) + k[1]*I
Ky = DiagonalMatrix(basis.kys) + k[2]*I
Kx_d = CuArray(Kx); Ky_d = CuArray(Ky); epc_d = CuArray(epc); muc_d = CuArray(muc);
if polarisation == TE
    LHS = Kx/epc*Kx + Ky/epc*Ky
    RHS = muc
    LHS_d = CuArray(LHS)
    RHS_d = CuArray(RHS)
    LHS_dd = Kx_d/epc_d*Kx_d + Ky_d/epc_d*Ky_d
    RHS_dd = muc_d
    label = L"H_z"
elseif polarisation == TM
    LHS = Kx/muc*Kx + Ky/muc*Ky
    RHS = epc
    label = L"E_z"
end
# Sometimes the generalised eigenvalue problem solver
# fails near Γ when the crystals are symmetric.
# In these cases, rewrite as a regular eigenvalue problem
@time freqs_squared, modes_data = eigen(RHS \ LHS);
@time freqs_squared, modes_data = eigen(Array(RHS_d \ LHS_d));
