using LinearAlgebra
using ProgressMeter

@enum Polarisation TE TM

function normalise(us; w=I)
    out = copy(us)
    for n in 1:size(out,2)
        out[:,n] /= sqrt(abs(dot(out[:,n],w*out[:,n])))
    end
    return out
end

function orthonormalise(us; w=I)
    out = copy(us)
    out[:,1] = normalise(out[:,1], w=w)
    # Gram-Schmidt orthogonalisation
    proj = zeros(size(out,1),size(out,1))
    for n in 2:size(out,2)
        proj += out[:,n-1]*out[:,n-1]' * w
        out[:,n] = out[:,n] - proj*out[:,n]
        out[:,n] = normalise(out[:,n], w=w)
    end
    return out
end

function solve(cg::ConvolvedGeometry, k, pol::Polarisation; ns=:, make_orthonormal=false)
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
    Λs, us = eigen(LHS, RHS)
    # Sort by increasing frequency
    ws = sqrt.(Λs)
    idx = sortperm(ws, by=real)
    ws = ws[idx][ns]
    us = us[:,idx][:,ns]
    if make_orthonormal
        us = orthonormalise(us, w=RHS)
    else
        us = normalise(us, w=RHS)
    end
    return ws, us
end
