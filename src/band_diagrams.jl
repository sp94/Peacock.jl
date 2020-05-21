"""
    sample_path(ks; [labels=[]], [dk=nothing])

Return a path through `ks` where `dk` is the maximum distance between points.

If `dk==nothing` then there will be approximately 10 points between
ks[1] and ks[2].
"""
function sample_path(ks; labels=[], dk=nothing)
    # Get labels of the corners
    if labels == []
        labels = ["" for x in ks]
    end
    @assert length(ks) == length(labels)
    # Replace 0 with [0,0]
    ks = [k == 0 ? [0,0] : k for k in ks]
    # Set default k-point density
    if dk == nothing
        dk = norm(ks[2]-ks[1])/10
    end
    # Sample between the corners such that spacing between points
    # is smaller than or equal to dk
    ks_sampled = [ks[1]]
    labels_sampled = [labels[1]]
    for i in 1:length(ks)-1
        k1 = ks[i]
        k2 = ks[i+1]
        d = norm(k2-k1)
        N = ceil(Int, d/dk)
        for n in 1:N
            k = k1 + n*(k2-k1)/N
            push!(ks_sampled, k)
            push!(labels_sampled, "")
        end
        push!(labels_sampled, labels[i+1])
    end
    return ks_sampled, labels_sampled
end


"""
    plot_band_diagram(my_solve::Function, ks, <keyword arguments>)

Plot the bands generated by `my_solve(k)` along `ks`.

Keyword arguments
- `dk=nothing`: maximum distance between points
- `labels=[]`: overwrite the labels for each `k` in `ks`
- `bands=(:)`: indices of the bands to plot
- `frequency_scale=1`: rescales the frequencies before plotting
- `color="k"`: color of the bands
- `markersize=nothing`: overwrite the size of each point
- `show_vlines=true`: plot vertical lines at each `k` in `ks`
"""
function plot_band_diagram(my_solve::Function, ks;
                    dk=nothing, labels=[], bands=:, frequency_scale=1,
                    color="k", markersize=nothing, show_vlines=true)
    # Add labels
    xs = cumsum([0; norm.(diff(ks))])
    xticks(xs, labels)
    if show_vlines
        for (x,label) in zip(xs,labels)
            if label != ""
                axvline(x, color="k", lw=1)
            end
        end
    end
    # Sample path
    ks, labels = sample_path(ks, labels=labels, dk=dk)
    xs = cumsum([0; norm.(diff(ks))])
    # Solve and plot points
    for (x,k) in zip(xs,ks)
        frequencies = my_solve(k)
        xs_k = [x for _ in 1:length(frequencies[bands])]
        ys = frequency_scale * frequencies[bands]
        PyPlot.plot(xs_k, real(ys), ".", color=color, alpha=0.5, markersize=markersize)
    end
    # Limits
    xlim(xs[1], xs[end])
    ylim(bottom=0)
end


"""
    plot_band_diagram(solver::Solver, ks, polarisation::Polarisation, <keyword arguments>)

Plot the bands generated by `solve(solver, k, polarisation)` along `ks`.
"""
function plot_band_diagram(solver::Solver, ks, polarisation::Polarisation;
            dk=0, labels=[], bands=1:10, frequency_scale=1, color="k", markersize=nothing)
    # Convert BrillouinZoneCoordinate to labelled positions in k space
    if labels == []
        labels = [hasproperty(x,:label) ? x.label : "" for x in ks]
    end
    ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]
    # Wrap all the variables into a single function of k that returns frequencies
    function my_solve(k)
        modes = solve(solver, k, polarisation)
        return [mode.frequency for mode in modes]
    end
    # Pass on to more general function
    plot_band_diagram(my_solve, ks; dk=dk, labels=labels, bands=bands, frequency_scale=frequency_scale, color=color, markersize=markersize)
end
