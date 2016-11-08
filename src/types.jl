import Base.zero
"""
Type to specify what bootstrap to run, and to hold the resulting bootstrapped values
"""
abstract AbstractBootstrapper{T}

type GenericBootstrapper{T} <: AbstractBootstrapper{T}
    nbootstrap::Int64
    func::Function
    values::Array{T,1}
end

type Bootstrapper{T<:Real} <: AbstractBootstrapper{T}
    nbootstrap::Int64
    func::Function
    values::Array{T,1}
    μ::T
    σ::T
end

function Bootstrapper{T}(::Type{T},nbootstrap::Int64,func::Function)
    #test the function
    a = func(rand(T,2))
    if length(a) > 1
        throw(ArgumentError("func should return a scalar value"))
    end
    Bootstrapper(nbootstrap,func,zeros(T,nbootstrap),zero(T),zero(T))
end

function summarize!{T<:Real}(B::AbstractBootstrapper{T})
    B.μ = mean(B.values)
    B.σ = std(B.values)
end

"""
Compute bootstrap statistics `nbootstrap` times using `func`

    function run_bootstrap{T}(nbootstrap::Int64,func::Function, X::Array{T,1},args...)
"""
function run_bootstrap{T}(nbootstrap::Int64,func::Function, X::Array{T,1},args...)
    B = Bootstrapper(T, nbootstrap, func)
    run_bootstrap!(B,X,args...)
    B
end

function run_bootstrap!{T}(B::AbstractBootstrapper{T}, X::Array{T,1},args...)
    n = length(X)
    idx = zeros(Int64,n)
    for i in 1:B.nbootstrap
        rand!(idx,1:n)
        B.values[i] = B.func(X[idx],args...)
    end
    summarize!(B)
end
