module StatUtils
using SpecialFunctions
using Random
using Statistics
using Optim
import StatsBase
using Distributions
import Distributions:pdf,logpdf, mean,var, fit_mle
include("types.jl")
include("generalized_poisson.jl")
export ZScore,Percentile, GeneralizedPoisson

"""
Compute the spearman rank correlation between r `x1` and `x2`. Also returns a two sided p-value for the signifiance of r using a permutation test
"""
function spearmanr(x1,x2,tail::Symbol=:right, RNG=MersenneTwister(rand(UInt32)))
    if !(tail in [:left, :right, :both])
        ArgumentError("tail must be one of :left, :right, :both")
    end
    cc = StatsBase.corspearman(x1,x2)
    nnl = 0
    nns = 0
    _x2 = copy(x2)
    N = 1000
    for i in 1:N
        _cc = StatsBase.corspearman(x1,shuffle!(RNG, _x2))
        if _cc > cc
            nnl +=1
        elseif _cc < cc
            nns +=1
        end
    end
    if tail == :both
        pv =  (nns + nnl)/N
    elseif tail == :left
        pv = nns/N
    else
        pv = nnl/N
    end
    cc, pv
end

"""
Use the False Discovery Rate (FDR) procedure in

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing JR Stat Soc B 57: 289–300.

to correct for the multiple comparisons represented by `pvalues` at level `0 < q < 1.0`.
"""
function fdr(pvalues::Array{Float64,1},q=0.05)
	sidx = sortperm(pvalues)
	m = length(sidx)
	vv = [1:m;]
	idx = findlast(q*vv./m .> pvalues[sidx])
	sidx[1:idx]
end

"""
Group the elements of `A` according to `grouping`, returning a dictionary of type `Dict{eltype(grouping), typeof(A)}`, where each entry (k,v) is a vector `v` with the elements of `A` corresponding to group `k`.

	function groupby{T1<:Any,T2<:Any}(A::AbstractArray{T1,1}, grouping::AbstractArray{T2,1})
"""
function groupby(A::AbstractArray{T1,1}, grouping::AbstractArray{T2,1}) where T1<:Any where T2<:Any
	groups = Dict{T2, Array{T1,1}}()
	for (a,g) in zip(A,grouping)
		if !(g in keys(groups))
			groups[g] = T1[]
		end
		push!(groups[g], a)
	end
	groups
end

function bootstrap_median(x1::AbstractArray{T,2}, n=10_000, RNG=MersenneTwister(rand(UInt32)),nx=size(x1,2)) where T<: Real
    nbins, ntrials = size(x1)
    μ = zeros(nbins)
    σ² = fill!(similar(μ), 0.0)
    nx = min(nx, ntrials)
    xs = zeros(size(x1,1), nx)
    for i in 1:n
        for k in 1:nx
            idx1 = rand(RNG, 1:ntrials)
            for j in 1:nbins
                xs[j,k] = x1[j,idx1]
            end
        end
        for j in 1:nbins
            m = median(xs[j,:])
            μ[j] += m
            σ²[j] += m*m
        end
    end
    μ ./= n
    σ² ./= n
    σ = sqrt.(σ² - μ.*μ)
    μ, σ
end

"""
Bootstrap linear regression of `y` on `x` by performing regression `n` times on random samples from `x` and `y`.
"""
function bootstrap_regression(x::AbstractVector{T},y::AbstractVector{T},n=1000;RNG=MersenneTwister(rand(UInt32))) where T <: Real
    nx = length(x)
    ym = Dict{T,T}()
    ys = Dict{T,T}()
    nn = Dict{T,Int64}()
    for i in 1:n
        _idx = rand(RNG, 1:nx, nx)
        _x = x[_idx]
        a,b = hcat(fill!(similar(_x), 1.0), _x)\(y[_idx])
        ye = a .+ b*_x
        for (_xe, _ye) in zip(_x, ye)
            ym[_xe] = get(ym, _xe, 0.0) + _ye
            ys[_xe] = get(ys, _xe, 0.0) + _ye*_ye
            nn[_xe]  = get(nn, _xe, 0) + 1
        end
    end
    xx = sort(collect(keys(ym)))
    μ = [ym[_xe]/nn[_xe] for _xe in xx]
    σ = sqrt.([ys[_xe]/nn[_xe] - (ym[_xe]/nn[_xe])^2 for _xe in xx])
    μ, σ, xx
end

"""
Regress `y` onto `x` by using the L₁ norm
"""
function robust_regression(x::AbstractVector{T}, y::AbstractVector{T},β0=rand(Float64,2)) where T <: Real
    xt = [fill!(similar(x), one(T)) x]
    func(β) = sum(abs, y .- xt*β)
    q = optimize(func, β0) 
    q.minimizer
end

"""
Compute the mean of each group as indicated by `grouping`.
"""
function Statistics.mean(x::AbstractArray{T}, grouping::AbstractVector{T2}) where T <: Real where T2 <: Integer
    ng = length(grouping)
    sx = size(x)
    dim = findfirst(sx.==ng)
    spre = sx[1:dim-1]
    spost = sx[dim+1:end]
    if dim == nothing
        error("`grouping` is not compatible with any dimension of `x`")
    end
    _grouping = compress(grouping)
    ngroups = maximum(_grouping)
    outsize = (sx[1:dim-1]..., ngroups, sx[dim+1:end]...)
    μ = fill(zero(T),outsize)
    counts = fill(0, outsize)
    for Ipost in CartesianIndices(spost)
        for i in axes(x,dim)
            for Ipre in CartesianIndices(spre)
                μ[Ipre,_grouping[i], Ipost] += x[Ipre, i, Ipost]
                counts[Ipre, _grouping[i], Ipost] += 1
            end
        end
    end
    μ ./= counts
    μ
end

"""
Compress the values in `x` by finding the smallest range that encompases all of `x` and repacing each value in `x` by its index into this range.
"""
function compress(x::AbstractArray{T}) where T <: Integer
    groups = unique(x)
    sort!(groups)
    y = fill!(similar(x), zero(T))
    for (i,_x) in enumerate(x)
        y[i] = findfirst(g->_x==g, groups)
    end
    y
end

end
