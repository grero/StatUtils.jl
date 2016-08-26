module StatUtils
import StatsBase
include("types.jl")
"""
Compute the spearman rank correlation between r `x1` and `x2`. Also returns a two sided p-value for the signifiance of r using a permutation test
"""
function spearmanr(x1,x2,tail::Symbol=:right)
    if !(tail in [:left, :right, :both])
        ArgumentError("tail must be one of :left, :right, :both")
    end
    cc = StatsBase.corspearman(x1,x2)
    nnl = 0
    nns = 0
    _x2 = copy(x2)
    N = 1000
    for i in 1:N
        _cc = StatsBase.corspearman(x1,shuffle!(_x2))
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
    pv
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
function groupby{T1<:Any,T2<:Any}(A::AbstractArray{T1,1}, grouping::AbstractArray{T2,1})
	groups = Dict{T2, Array{T1,1}}()
	for (a,g) in zip(A,grouping)
		if !(g in keys(groups))
			groups[g] = T1[]
		end
		push!(groups[g], a)
	end
	groups
end

end

