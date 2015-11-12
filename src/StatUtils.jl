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

end

