using LinearAlgebra
using ProgressMeter

@enum Polarisation TE TM

function orthogonalise(us; w=I)
    # Gram-Schmidt orthogonalisation
    proj = zeros(size(us,1),size(us,1))
    for n in 2:size(us,2)
        proj += us[:,n-1]*us[:,n-1]' * w
        us[:,n] = us[:,n] - proj*us[:,n]
        us[:,n] /= sqrt(abs(dot(us[:,n],w*us[:,n])))
    end
end

function solve(cg::ConvolvedGeometry, k, pol::Polarisation)
	# Eigenvalue problem
    epc, muc = cg.epc, cg.muc
    kx = cg.kx + k[1]*I
    ky = cg.ky + k[2]*I
    if pol == TM
        LHS = kx/muc*kx + ky/muc*ky
        RHS = epc
    elseif pol == TE
        LHS = kx/epc*kx + ky/epc*ky
        RHS = muc
    end
    vals, vecs = eigen(LHS, RHS)
    vecs = orthogonalise(vecs, w=RHS)
    ws = sqrt.(vals)
    us = vecs
    # Sort by increasing frequency
    idx = sortperm(ws, by=real)
    ws = ws[idx]
    us = us[:,idx]
    return ws, us
end
