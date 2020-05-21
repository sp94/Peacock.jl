using LinearAlgebra


"""
    overlaps(a::HilbertSpace, b::HilbertSpace)

Calculate the overlaps between the basis vectors of each Hilbert space.
"""
function overlaps(a::HilbertSpace, b::HilbertSpace)
    @assert a.weighting == b.weighting
    LHS = a.data
    RHS = b.weighting * b.data
    return [dot(LHS[:,i],RHS[:,j]) for i in 1:size(a.data,2), j in 1:size(b.data,2)]
end


"""
    unitary_approx(M::AbstractArray)

Calculate the best unitary approxmation of `M` using singular value
decomposition.
"""
function unitary_approx(M::AbstractArray)
    F = svd(M)
    return F.U * F.Vt
end


"""
    unitary_overlaps(a::HilbertSpace, b::HilbertSpace)

Shortcut for unitary_approx(overlaps(a, b))
"""
function unitary_overlaps(a::HilbertSpace, b::HilbertSpace)
    unitary_overlaps = overlaps(a, b)
    return unitary_approx(unitary_overlaps)
end


"""
    wilson_matrix(spaces::Array{HilbertSpace,1}; closed::Bool=true)

Calculate the Wilson loop matrix through the Hilbert spaces.

The `closed` keyword 
"""
function wilson_matrix(spaces::Array{HilbertSpace,1}; closed::Bool=true)
    @assert length(spaces) > 1
    # Assume all spaces have the same plane wave expansion
    basis = spaces[1].basis
    for space in spaces
        @assert space.basis == basis
    end
    # However, k0 is still changing...
    if closed
        # In a closed Wilson loop...
        @assert spaces[1] == spaces[end]
    else
        # If not closed, then the loop begins and finishes at
        # the same k0 but in different Brillouin zones
        delta_k0 = spaces[end].k0 - spaces[1].k0
        b1 = spaces[1].basis.b1
        b2 = spaces[1].basis.b2
        B = [b1 b2]
        dp, dq = B \ delta_k0
        # Round to integers if approximately correct
        if isapprox(round(dp), dp, atol=1e-6) && isapprox(round(dq), dq, atol=1e-6)
            dp = Int(round(dp))
            dq = Int(round(dq))
        else
            @show dp, dq
            throw("First and last spaces do not... (TODO)")
        end
        # Duplicate spaces[1] and shift its k0 so that it matches spaces[end]
        space_1_at_end = shift_k0(spaces[1], dp, dq)
        spaces = [spaces; space_1_at_end]
    end
    W = I
    for (a,b) in zip(spaces, spaces[2:end])
        W = W * unitary_overlaps(a, b)
    end
    @assert W * W' ≈ I
    return W
end


"""
    wilson_eigvals(spaces::AbstractArray{HilbertSpace,1}; closed=true)

Return the eigenvalues of the Wilson matrix through a Hilbert space,
sorted by phase angle.
"""
function wilson_eigvals(spaces::AbstractArray{HilbertSpace,1}; closed=true)
    W = wilson_matrix(spaces, closed=closed)
    vals = eigvals(W)
    return sort(vals, by=angle)
end


"""
    wilson_eigvals(spaces::AbstractArray{HilbertSpace,1}; closed=true)

Return the eigenvalues and eigenvectors of the Wilson matrix through a Hilbert
space, sorted by the phase angle of the eigenvalues.
"""
function wilson_eigen(spaces::AbstractArray{HilbertSpace,1}; closed=true)
    W = wilson_matrix(spaces, closed=closed)
    vals, vecs = eigen(W)
    idx = sortperm(vals, by=angle)
    vals = vals[idx]
    vecs = vecs[:,idx]
    vecs = unitary_approx(vecs) # or other orthonormalisation? 
    return vals, vecs
end


"""
    wilson_gauge(spaces::AbstractArray{HilbertSpace,1}; closed=true)

Return the eigenvalues, eigenvectors, and gauge of the Wilson loop through
the Hilbert space, sorted by the phase angle of the eigenvalues.
"""
function wilson_gauge(spaces::AbstractArray{HilbertSpace,1}; closed=true)
    vals, vecs = wilson_eigen(spaces, closed=closed)
    gauge = copy(spaces)
    gauge[1] = HilbertSpace(gauge[1].k0, gauge[1].data*vecs, gauge[1].weighting, gauge[1].basis)
    for i in 1:length(gauge)-1
        a = gauge[i]
        b = gauge[i+1]
        mixing = unitary_overlaps(b, a)
        gauge[i+1] = HilbertSpace(b.k0, b.data*mixing, b.weighting, b.basis)
    end
    if closed
        @assert overlaps(gauge[end], gauge[1]) ≈ diagm(0=>vals)
    end
    return vals, vecs, gauge
end


"""
    plot_wilson_loop_winding(solver::Solver, ks, polarisation, bands::AbstractVector{<:Int}, <keyword arguments>)

Perform a series of Wilson loops along `ks`, and plot the spectra on a band diagram.

Keyword arguments
- `dk=nothing`: maximum distance between points
- `labels=[]`: overwrite the labels for each `k` in `ks`
- `delta_brillouin_zone=BrillouinZoneCoordinate(0,1)`: each Wilson loop starts at and finishes in at the same `k` in different Brillouin zones
"""
function plot_wilson_loop_winding(solver::Solver, ks, polarisation, bands::AbstractVector{<:Int};
                dk=0, labels=[], delta_brillouin_zone=BrillouinZoneCoordinate(0,1))
    # Convert BrillouinZoneCoordinate to labelled positions in k space
    if labels == []
        labels = [hasproperty(x,:label) ? x.label : "" for x in ks]
    end
    ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]
    # Function to return the Wilson loop eigenvalues
    function my_solve(k)
        spaces = HilbertSpace[]
        delta_k = get_k(delta_brillouin_zone, solver.basis)
        ks_inner = [k, k+delta_k]
        ks_inner, _ = sample_path(ks_inner, dk=dk)
        for k_inner in ks_inner
            modes = solve(solver, k_inner, polarisation)
            space = HilbertSpace(modes[bands])
            push!(spaces, space)
        end
        vals = wilson_eigvals(spaces, closed=false)
        angles = angle.(vals)
        angles = [angles.-2pi angles angles.+2pi]
        return angles
    end
    plot_band_diagram(my_solve, ks, dk=dk, labels=labels, show_vlines=false)
    ylim(-pi,pi)
    yticks((-1:1)*pi, labels=["-π","0","+π"])
    ylabel("Wilson spectrum")
end
