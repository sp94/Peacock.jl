using LinearAlgebra
using ProgressMeter

@enum Polarisation TE TM

function orthogonalise(us; w=I)
    out = copy(us)
    out[:,1] /= sqrt(abs(dot(out[:,1],w*out[:,1])))
    # Gram-Schmidt orthogonalisation
    proj = zeros(size(out,1),size(out,1))
    for n in 2:size(out,2)
        proj += out[:,n-1]*out[:,n-1]' * w
        out[:,n] = out[:,n] - proj*out[:,n]
        out[:,n] /= sqrt(abs(dot(out[:,n],w*out[:,n])))
    end
    return out
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
    # Sort by increasing frequency
    ws = sqrt.(vals)
    idx = sortperm(ws, by=real)
    ws = ws[idx]
    vecs = vecs[:,idx]
    # Ensure orthogonal, as we have a Hermitian eigenvalue problem
    us = orthogonalise(vecs, w=RHS)
    return ws, us
end
