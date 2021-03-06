import Base.zero, Base.rand, Base.convert, Base.abs, Base.typemax, Base.typemin
"""
Type to specify what bootstrap to run, and to hold the resulting bootstrapped values
"""
abstract type AbstractBootstrapper{T} end

mutable struct GenericBootstrapper{T} <: AbstractBootstrapper{T}
    nbootstrap::Int64
    func::Function
    values::Array{T,1}
end

mutable struct Bootstrapper{T<:Real} <: AbstractBootstrapper{T}
    nbootstrap::Int64
    func::Function
    values::Array{T,1}
    μ::T
    σ::T
end

function Bootstrapper(::Type{T},nbootstrap::Int64,func::Function) where T <: Real
    #test the function
    a = func(rand(T,2))
    if length(a) > 1
        throw(ArgumentError("func should return a scalar value"))
    end
    Bootstrapper(nbootstrap,func,zeros(T,nbootstrap),zero(T),zero(T))
end

function summarize!(B::AbstractBootstrapper{T}) where T<:Real
    B.μ = mean(B.values)
    B.σ = std(B.values)
end

"""
Compute bootstrap statistics `nbootstrap` times using `func`

    function run_bootstrap{T}(nbootstrap::Int64,func::Function, X::Array{T,1},args...)
"""
function run_bootstrap(nbootstrap::Int64,func::Function, X::Array{T,1},args...;RNG=MersenneTwister(rand(UInt32))) where T <: Any
    B = Bootstrapper(T, nbootstrap, func)
    run_bootstrap!(B,X,args...;RNG=RNG)
    B
end

function run_bootstrap!(B::AbstractBootstrapper{T}, X::Array{T,1},args...;RNG=MersenneTwister(rand(UInt32))) where T <: Any
    n = length(X)
    idx = zeros(Int64,n)
    for i in 1:B.nbootstrap
        rand!(RNG, idx,1:n)
        B.values[i] = B.func(X[idx],args...)
    end
    summarize!(B)
end

struct ZScore
	v::Float64
end

zero(::Type{ZScore}) = ZScore(zero(Float64))
rand(::Type{ZScore}) = randn()
convert(::Type{ZScore}, X::Float64) = ZScore(X)
convert(::Type{Float64}, X::ZScore) = X.v
abs(z::ZScore) = abs(z.v)

struct Percentile
    v::Float64
    Percentile(v) = 0.0 < v < 100.0 ? new(v) : throw(ArgumentError("Percentiles must be between 0.0 and 100.0"))
end

Base.typemax(::Type{Percentile}) = 100.0
Base.typemin(::Type{Percentile}) = 0.0
