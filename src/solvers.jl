using LinearAlgebra
using ProgressMeter

@enum Polarisation TE TM

function solve(cg::ConvolvedGeometry, k, pol::Polarisation)
	# Eigenvalue problem
    epc, muc = cg.epc, cg.muc
    kx = cg.kx + k[1]*I
    ky = cg.ky + k[2]*I
    if pol == TM
        eigvals, eigvecs = eigen(kx/muc*kx+ky/muc*ky, epc)
    elseif pol == TE
        eigvals, eigvecs = eigen(kx/epc*kx+ky/epc*ky, muc)
    end
    ws = sqrt.(eigvals)
    us = eigvecs
    # Normalisation 
    for j in 1:size(us,2)
        if pol == TM
            # TM modes are orthogonal under <ui|ep|uj> = δ_ij
    	    us[:,j] /= sqrt(abs(dot(us[:,j],epc*us[:,j])))
        elseif pol == TE
            # TE modes are orthogonal under <ui|mu|uj> = δ_ij
            us[:,j] /= sqrt(abs(dot(us[:,j],muc*us[:,j])))
        end
    end
    # Sort by increasing frequency
    idx = sortperm(ws, by=real)
    ws = ws[idx]
    us = us[:,idx]
    return ws, us
end
