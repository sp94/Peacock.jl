function sample_path(k_path; dk=0)
	# Replace 0 with [0,0]
	k_path = Array{Float64,1}[k == 0 ? [0,0] : k for k in k_path]
	# Set default k-point density
	if dk == 0
		dk = norm(k_path[2]-k_path[1])/10
	end
	out = [k_path[1]]
	for (k1,k2) in zip(k_path, k_path[2:end])
		# Sample the region between k1 and k2 such that
		# the spacing between points is smaller than or equal to dk
		d = norm(k2-k1)
		N = ceil(d/dk)
		for n in 1:N
			k = k1 + n*(k2-k1)/N
			push!(out, k)
		end
	end
	return out
end