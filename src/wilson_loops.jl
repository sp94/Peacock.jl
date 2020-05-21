using LinearAlgebra


"""
To do
"""
function overlaps(a::HilbertSpace, b::HilbertSpace)
    @assert a.weighting == b.weighting
    LHS = a.data
    RHS = b.weighting * b.data
    return [dot(LHS[:,i],RHS[:,j]) for i in 1:size(a.data,2), j in 1:size(b.data,2)]
end


"""
To do
"""
function unitary_approx(M::AbstractArray)
    F = svd(M)
    return F.U * F.Vt
end


"""
To do
"""
function unitary_overlaps(a::HilbertSpace, b::HilbertSpace)
    unitary_overlaps = overlaps(a, b)
    return unitary_approx(unitary_overlaps)
end

##### move into tests?
# function assert_orthonormal(a::Eigenvectors)
#     x = overlaps(a, a)
#     for p in 1:size(x,1), q in 1:size(x,2)
#         @assert isapprox(x[p,q], p==q, atol=1e-6)
#     end
# end


"""
To do
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
To do
"""
function wilson_eigvals(spaces::AbstractArray{HilbertSpace,1}; closed=true)
    W = wilson_matrix(spaces, closed=closed)
    vals = eigvals(W)
    return sort(vals, by=angle)
end


"""
To do
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
To do
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
Provide a wrapper around plot_band_diagram that is focused on plotting Wilson loops
"""
function plot_wilson_loop_winding(solver::Solver, ks, polarisation, bands::AbstractVector{<:Int}; dk=0, labels=[],
            loop_direction=BrillouinZoneCoordinate(0,1), ts=range(0,stop=1,length=21))
    # Convert BrillouinZoneCoordinate to labelled positions in k space
    if labels == []
        labels = [hasproperty(x,:label) ? x.label : "" for x in ks]
    end
    ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]
    # Function to return the Wilson loop eigenvalues
    function my_solve(k)
        spaces = HilbertSpace[]
        for t in ts
            dk = t * get_k(loop_direction, solver.basis)
            freqs, modes = solve(solver, k+dk, polarisation)
            space = HilbertSpace(modes[bands])
            push!(spaces, space)
        end
        vals = wilson_eigvals(spaces, closed=false)
        angles = angle.(vals)
        angles = [angles.-2pi angles angles.+2pi]
        return angles
    end
    @show dk
    plot_band_diagram(my_solve, ks, dk=dk, labels=labels, show_vlines=false)
    ylim(-pi,pi)
    yticks((-1:1)*pi, labels=["-π","0","+π"])
    ylabel("Wilson spectrum")
end
