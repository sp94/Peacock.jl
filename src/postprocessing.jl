"""
Change the basis of the plane wave expansion.
This should correspond to a symmetry of the system,
eg shift by a reciprocal lattice vector or a rotation
"""
function change_basis(cg, us, pq_map)
    us_ = zeros(ComplexF64, size(us))
    for row in 1:size(us,1)
        # Identify the corresponding (p_,q_) indices
        # in the new plane wave basis
        p = cg.ps[row]
        q = cg.ps[row]
        p_, q_ = pq_map(p, q)
        # Copy only if the plane wave still exists in the new basis
        row_ = findfirst(x->x==(p_,q_), collect(zip(cg.ps,cg.qs)))
        if !isnothing(row_)
            us_[row_,:] = us[row,:]
        end
    end
    return us_
end


"""
Shift the plane wave expansion by a reciprocal lattice vector, p*b1+q*b2.
"""
function shift(cg, us, dp::Int, dq::Int)
    pq_map(p,q) = (p-dp, q-dq)
    return change_basis(cg, us, pq_map)
end


"""
Convert plane wave expansion to real space field.

    get_field(cg, u; k0=[0,0], t1s=-0.5:0.01:0.5, t2s=-0.5:0.01:0.5)

Leaving k0=[0,0] will return the cell-periodic part of the Bloch field only.
Otherwise, k0 must be the wavevector for which u was calculated.

The field will have dimensions length(t1s) x length(t2s).
"""
function get_field(cg, u; k0=[0,0], t1s=-0.5:0.01:0.5, t2s=-0.5:0.01:0.5)
    # Generate real space coordinates
    xs = [t1*cg.a1[1]+t2*cg.a2[1] for t1 in t1s, t2 in t2s]
    ys = [t1*cg.a1[2]+t2*cg.a2[2] for t1 in t1s, t2 in t2s]
    out = zeros(ComplexF64, length(t1s), length(t2s))
    # Inverse Fourier transform
    for (p,q,u_pq) in zip(cg.ps, cg.qs, u)
        k = k0 + p*cg.b1 + q*cg.b2
        for c in CartesianIndices(out)
            bloch_phase = k[1]*xs[c] + k[2]*ys[c]
            out[c] += u_pq * exp(1im*bloch_phase)
        end
    end
    return out
end
