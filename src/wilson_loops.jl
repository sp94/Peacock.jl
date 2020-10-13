using LinearAlgebra


"""
    overlaps(a::Eigenspace, b::Eigenspace)

Calculate the overlaps between the basis vectors of each eigenspace.
"""
function overlaps(a::Eigenspace, b::Eigenspace)
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
    unitary_overlaps(a::Eigenspace, b::Eigenspace)

Shortcut for unitary_approx(overlaps(a, b))
"""
function unitary_overlaps(a::Eigenspace, b::Eigenspace)
    unitary_overlaps = overlaps(a, b)
    return unitary_approx(unitary_overlaps)
end


"""
    wilson_matrix(spaces::Array{Eigenspace,1}; k_map=nothing)

Calculate the Wilson loop matrix through the eigenspaces.

The `k_map` function should describe the transformation applied to `spaces[1]`
to bring it to the end of the loop, ie, `k_map(spaces[1]) == spaces[end].k0`
This could be a rotation, translation, reflection, etc, in k-space.

If the `k_map` keyword is not set, it will be assumed that the first
and last k-points are connected by a translation in k-space.
"""
function wilson_matrix(spaces::Array{Eigenspace,1}; k_map=nothing)
    @assert length(spaces) > 1
    # Assume all spaces have the same plane wave expansion
    basis = spaces[1].basis
    for space in spaces
        @assert space.basis == basis
    end
    # However, k0 is still changing...
    if k_map == nothing
        # If the user hasn't specified a k_map, then
        # assume Wilson loop is closed by a translation in k-space
        _k_map(k0) = k0 + spaces[end].k0 - spaces[1].k0
    else
        # Otherwise use their k_map
        _k_map = k_map
    end
    @assert _k_map(spaces[1].k0) ≈ spaces[end].k0
    space_1_at_end = transform(spaces[1], _k_map)
    spaces = [spaces; space_1_at_end]
    W = I
    for (a,b) in zip(spaces, spaces[2:end])
        W = W * unitary_overlaps(a, b)
    end
    @assert W * W' ≈ I
    return W
end


"""
    wilson_eigvals(spaces::AbstractArray{Eigenspace,1}; k_map=nothing)

Return the eigenvalues of [`wilson_matrix`](@ref), sorted by phase angle.
"""
function wilson_eigvals(spaces::AbstractArray{Eigenspace,1}; k_map=nothing)
    W = wilson_matrix(spaces, k_map=k_map)
    vals = eigvals(W)
    return sort(vals, by=angle)
end


"""
    wilson_eigvals(spaces::AbstractArray{Eigenspace,1}; k_map=nothing)

Return the eigenvalues and eigenvectors of [`wilson_matrix`](@ref),
sorted by the phase angle of the eigenvalues.
"""
function wilson_eigen(spaces::AbstractArray{Eigenspace,1}; k_map=nothing)
    W = wilson_matrix(spaces, k_map=k_map)
    vals, vecs = eigen(W)
    idx = sortperm(vals, by=angle)
    vals = vals[idx]
    vecs = vecs[:,idx]
    vecs = unitary_approx(vecs) # or other orthonormalisation? 
    return vals, vecs
end


"""
    wilson_gauge(spaces::AbstractArray{Eigenspace,1}; k_map=nothing)

Return the eigenvalues, eigenvectors, and gauge of the Wilson loop through
the eigenspace, sorted by the phase angle of the eigenvalues.
"""
function wilson_gauge(spaces::AbstractArray{Eigenspace,1}; k_map=nothing)
    vals, vecs = wilson_eigen(spaces, k_map=k_map)
    gauge = copy(spaces)
    gauge[1] = Eigenspace(gauge[1].k0, gauge[1].data*vecs, gauge[1].weighting, gauge[1].basis, gauge[1].data_label)
    for i in 1:length(gauge)-1
        a = gauge[i]
        b = gauge[i+1]
        mixing = unitary_overlaps(b, a)
        gauge[i+1] = Eigenspace(b.k0, b.data*mixing, b.weighting, b.basis, b.data_label)
    end
    return vals, vecs, gauge
end


"""
    plot_wilson_loop_winding(solver::Solver, ks, polarisation, bands::AbstractVector{<:Int}, <keyword arguments>)

Perform a series of Wilson loops along `ks`, and plot the spectra on a band diagram.

Keyword arguments
- `dk_outer=nothing`: maximum distance between each loop (resolution of the scan)
- `dk_inner=nothing`: maximum distance between points along a loop (resolution of the loop)
- `delta_brillouin_zone=BrillouinZoneCoordinate(0,1)`: each Wilson loop starts at and finishes in at the same `k` in different Brillouin zones
- `labels=[]`: overwrite the labels for each `k` in `ks`
- `markersize=nothing`: overwrite the size of each point
"""
function plot_wilson_loop_winding(solver::Solver, ks, polarisation, bands::AbstractVector{<:Int};
                dk_outer=nothing, dk_inner=nothing, delta_brillouin_zone=BrillouinZoneCoordinate(0,1),
                labels=[], markersize=nothing)
    # Convert BrillouinZoneCoordinate to labelled positions in k space
    if labels == []
        labels = [hasproperty(x,:label) ? x.label : "" for x in ks]
    end
    ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]
    # Function to return the Wilson inner eigenvalues
    function my_solve(k)
        spaces = Eigenspace[]
        delta_k = get_k(delta_brillouin_zone, solver.basis)
        ks_inner = [k, k+delta_k]
        ks_inner, _ = sample_path(ks_inner, dk=dk_inner)
        for k_inner in ks_inner
            modes = solve(solver, k_inner, polarisation)
            space = Eigenspace(modes[bands])
            push!(spaces, space)
        end
        vals = wilson_eigvals(spaces)
        angles = angle.(vals)
        angles = [angles.-2pi angles angles.+2pi]
        return angles
    end
    plot_band_diagram(my_solve, ks, dk=dk_outer, labels=labels, show_vlines=false, markersize=markersize)
    ylim(-pi,pi)
    yticks((-1:1)*pi, labels=["-π","0","+π"])
    ylabel("Wilson spectrum")
end
