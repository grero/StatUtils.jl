struct GeneralizedPoisson <: Distributions.DiscreteUnivariateDistribution
    θ::Float64
    λ::Float64
end

function Distributions.pdf(G::GeneralizedPoisson, x::Int64)
    G.θ*(G.θ+G.λ*x)^(x-1)*exp(-G.θ - G.λ*x)/factorial(x)
end

function Distributions.logpdf(G::GeneralizedPoisson, x::Int64)
    a = log(G.θ)
    a += (x-1)*log(G.θ+G.λ*x)
    a -= G.θ
    a -= G.λ*x
    a -= lgamma(x+1)
    a
end

Distributions.mean(G::GeneralizedPoisson) = G.θ/(1-G.λ)
Distributions.var(G::GeneralizedPoisson) = G.θ/(1-G.λ)^3
