using LinearAlgebra
using ProgressMeter

@enum Polarisation TE TM

struct Solver
    cg::ConvolvedGeometry
    pol::Polarisation
end

function get_equation(solver::Solver, k::AbstractArray{<:Real,1})
    epc, muc = solver.cg.epc, solver.cg.muc
    kx = solver.cg.kx + k[1]*I
    ky = solver.cg.ky + k[2]*I
    if solver.pol == TM
        LHS = kx/muc*kx+ky/muc*ky
        RHS = epc
    elseif solver.pol == TE
        LHS = kx/epc*kx+ky/epc*ky
        RHS = muc
    end
    return LHS, RHS
end

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


struct Eigenvectors
    data::AbstractMatrix{ComplexF64}
    weighting::Union{UniformScaling,AbstractMatrix}
end


function solve(solver::Solver, k::AbstractVector{<:Real}; ns=:, make_orthonormal=false)
    # Solve generalised eigenvalue problem
    LHS, RHS = get_equation(solver, k)
    ws_squared, us_data = try
        eigen(LHS, RHS)
    catch
        # Sometimes the generalised eigenvalue problem solver
        # fails near Γ when the crystals are symmetric.
        # In these cases, rewrite as a regular eigenvalue problem
        eigen(RHS \ LHS)
    end
    # Sort by increasing frequency
    ws = sqrt.(ws_squared)
    idx = sortperm(ws, by=real)
    ws = ws[idx][ns]
    us_data = us_data[:,idx][:,ns]
    # Eigenvectors are weighted by the RHS of the generalised eigenvalue problem
    weighting = RHS
    if make_orthonormal
        us_data = orthonormalise(us_data, w=weighting)
    else
        us_data = normalise(us_data, w=weighting)
    end
    us = Eigenvectors(us_data, RHS)
    return ws, us
end
