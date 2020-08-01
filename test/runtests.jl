using Test
using Peacock
using LinearAlgebra

using Peacock.Zoo
using Parameters

# Test sample_path
begin
    # Generate a random path
    ks_orig = [[0.0,0.0]]
    for n in 1:3
        push!(ks_orig, ks_orig[end]+rand(2))
    end
    labels_orig = ["k$(n)" for n in 1:length(ks_orig)]
    # Sample it
    dk = 0.1
    ks, labels = Peacock.sample_path(ks_orig, labels=labels_orig, dk=dk)
    # Check same number of ks and labels
    @test length(ks) == length(labels)
    # Check k-point spacing is correct
    for (k1,k2) in zip(ks,ks[2:end])
        @test norm(k2-k1) <= dk
    end
    # Check labels match expected k positions
    for (k_orig,label_orig) in zip(ks_orig,labels_orig)
        n = findfirst([k==k_orig for k in ks])
        @test labels[n] == label_orig
    end
end


# Test convmat
begin
    # Create a test basis of plane waves
    b1, b2 = [0.3,0.2], [0.5,-0.1]
    ps, qs = Int[], Int[]
    for p in -1:2, q in -2:3
        push!(ps, p)
        push!(qs, q)
    end
    basis = Peacock.PlaneWaveBasis(b1, b2, ps, qs)
    # Convolution matrix of homogeneous array should be diagonal
    x = 1 + rand(ComplexF64)
    mat = fill(x, 10, 20)
    cmat = Peacock.convmat(mat, basis)
    @test all(diag(cmat) .≈ x)
end


# Test homogeneous systems with known solutions
begin
    # Create homogeneous geometry
    ep = 1 + rand()
    mu = 1 + rand()
    a1 = [1,0]
    a2 = [0,1]
    d1 = d2 = 0.01
    geometry = Geometry((x,y)->ep, (x,y)->mu, a1, a2, d1, d1)
    # Check for increasing Fourier space cutoffs
    for cutoff in [1,3,5]
        # test with norm(b1) == norm(b2)
        solver = Solver(geometry, cutoff)
        for polarisation in [TE,TM]
            modes = solve(solver, [1,0], polarisation)
            @test isapprox(modes[1].frequency, 1/sqrt(ep*mu))
        end
        # test with norm(b1) != norm(b2)
        solver = Solver(geometry, cutoff, cutoff+2)
        for polarisation in [TE,TM]
            modes = solve(solver, [1,0], polarisation)
            @test isapprox(modes[1].frequency, 1/sqrt(ep*mu))
        end
    end
end


# Test conversion between real space and reciprocal space lattice vectors
a1 = rand(2)
a2 = rand(2)
b1, b2 = Peacock.as_to_bs(a1,a2)
@test isapprox(dot(a1,b1), 2pi)
@test isapprox(dot(a2,b2), 2pi)
@test isapprox(dot(a1,b2), 0, atol=1e-6)
@test isapprox(dot(a2,b1), 0, atol=1e-6)
a1_, a2_ = Peacock.bs_to_as(b1, b2)
@test isapprox(a1, a1_)
@test isapprox(a2, a2_)


# Test normalisation
data = rand(ComplexF64, 5, 5)
weighting = I + 0.5rand(ComplexF64, 5, 5)
weighting = weighting + weighting' # random positive-definite Hermitian weighting
data_ = Peacock.normalise(data, weighting=weighting)
for n in 1:size(data_,2)
    weighted_norm = sqrt(abs(dot(data_[:,n], weighting*data_[:,n])))
    @test isapprox(weighted_norm, 1)
end

# Test orthonormalisation
data_ = Peacock.orthonormalise(data, weighting=weighting)
for i in 1:size(data_,2), j in 1:size(data_,2)
    weighted_overlap = sqrt(abs(dot(data_[:,i], weighting*data_[:,j])))
    @test isapprox(weighted_overlap, float(i==j), atol=1e-6)
end


# Test that the Wilson loop gauge and Wilson loop eigenvalues agree
begin
    # Load solver from the Zoo
    @unpack solver, polarisation = make_wu_topo(11)
    # Small Wilson loop of bands 1-3
    k0 = [0.3,0.2] # arbitrary
    ks = [k0, k0+[0.1,0], k0+[0.1,0.1], k0+[0,0.1]]
    spaces = HilbertSpace[]
    for k in ks
        modes = solve(solver, k, polarisation)
        space = HilbertSpace(modes[1:3])
        push!(spaces, space)
    end
    push!(spaces, spaces[1])
    vals, vecs, gauge = Peacock.wilson_gauge(spaces)
    @test Peacock.overlaps(gauge[end], gauge[1]) ≈ diagm(0=>vals)
end


# Test Wilson loops of fragile topological crystal from Blanco de Paz et al 2019
begin

    # Load solver from the Zoo
    @unpack solver, polarisation = make_dePaz_frag(11)

    # Create a Wilson loop of bands 2&3 from k0 to k0+1
    function test_wilson_eigvals(solver, x::BrillouinZoneCoordinate, polarisation)
        k0 = Peacock.get_k(x, solver.basis)
        ts = range(0, stop=1, length=20)
        ks = [k0 + t*solver.basis.b2 for t in ts]
        spaces = HilbertSpace[]
        for k in ks
            modes = solve(solver, k, polarisation)
            space = HilbertSpace(modes[2:3])
            push!(spaces, space)
        end
        return Peacock.wilson_eigvals(spaces, closed=false)
    end

    # Wilson loop of bands 2&3 from Γ to Γ+b1
    vals = test_wilson_eigvals(solver, BrillouinZoneCoordinate(0.0,0.0), polarisation)
    @show vals
    # Expected eigenvalues = [-1, -1]
    percentage_error = angle.(vals / -1) / 2pi * 100
    @test all(abs.(percentage_error) .< 1)

    # Wilson loop of bands 2&3 from M to M+b1
    vals = test_wilson_eigvals(solver, BrillouinZoneCoordinate(0.5,0.0), polarisation)
    @show vals
    # Expected eigenvalues = [+1, +1]
    percentage_error = angle.(vals / +1) / 2pi * 100
    @test all(abs.(percentage_error) .< 1)

end
