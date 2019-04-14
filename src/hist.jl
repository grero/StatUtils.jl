"""
Compute a histogram of the angles `θ` using equally sized bins spanning [-π,π]
"""
function angular_histogram(θ::AbstractVector{T}, Δθ::Real) where T <: Real
    nbins = round(Int64,2π/Δθ)
    bins = collect(range(0.0, step=Δθ,length=div(nbins,2)))
    append!(bins, range(-pi, step=Δθ, length=div(nbins,2)))
    nbins = length(bins)
    sort!(bins)
    counts = fill(0, nbins)
    binidx = fill(0, length(θ))
    for i in 1:nbins
        idx = bins[i] - Δθ/2 .<= θ .< bins[i] + Δθ/2
        counts[i] = sum(idx)
        binidx[idx] .= i
    end
    counts, binidx, bins
end
