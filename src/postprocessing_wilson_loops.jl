using LinearAlgebra

function overlaps(a::Eigenvectors, b::Eigenvectors)
    @assert a.weighting == b.weighting
    LHS = a.data
    RHS = b.weighting * b.data
    return [dot(LHS[:,i],RHS[:,j]) for i in 1:size(a.data,2), j in 1:size(b.data,2)]
end

function unitary_approx(M::AbstractArray)
    F = svd(M)
    return F.U * F.Vt
end

function unitary_overlaps(a, b; w=I)
    unitary_overlaps = overlaps(a, b)
    return unitary_approx(unitary_overlaps)
end

function assert_orthonormal(a::Eigenvectors)
    x = overlaps(a, a)
    for p in 1:size(x,1), q in 1:size(x,2)
        @assert isapprox(x[p,q], p==q, atol=1e-6)
    end
end


# Wilson loop matrix
function wilson_matrix(us_list; ns=:, closed=true)
    # Could add some 'tolerance' argument to sample the path
    # until the projectors are smooth
    # Failure could indicate that the user hasn't provided an
    # isolated set of bands
    @assert length(us_list) > 1
    if closed
        @assert us_list[1] == us_list[end]
    end
    for a in us_list
        assert_orthonormal(a)
    end
    W = I
    for (a,b) in zip(us_list,us_list[2:end])
        a = Eigenvectors(a.data[:,ns], a.weighting)
        b = Eigenvectors(b.data[:,ns], b.weighting)
        W = W * unitary_overlaps(a, b)
    end
    @show W
    @show -1im*log(W)
    @assert W * W' ≈ I
    return W
end

function wilson_eigvals(us_list; ns=:, closed=true)
    W = wilson_matrix(us_list, ns=ns, closed=closed)
    vals = eigvals(W)
    return sort(vals, by=angle)
end

function wilson_eigen(us_list; ns=:, closed=true)
    W = wilson_matrix(us_list, ns=ns, closed=closed)
    vals, vecs = eigen(W)
    idx = sortperm(vals, by=angle)
    vals = vals[idx]
    vecs = vecs[:,idx]
    vecs = unitary_approx(vecs) # or other orthonormalisation? 
    return vals, vecs
end

function wilson_gauge(us_list; ns=:,closed=true)
    vals, vecs = wilson_eigen(us_list, ns=ns, closed=closed)
    gauge = copy(us_list)
    gauge[1] = Eigenvectors(gauge[1].data*vecs, gauge[1].weighting)
    for i in 1:length(gauge)-1
        a = gauge[i]
        b = gauge[i+1]
        mixing = unitary_overlaps(b, a)
        gauge[i+1] = Eigenvectors(b.data*mixing, b.weighting)
    end
    @assert overlaps(gauge[end], gauge[1]) ≈ diagm(0=>vals)
    return vals, vecs, gauge
end
