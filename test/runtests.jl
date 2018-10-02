using Distributions
using StatUtils
using Test
using Random

function test_bootstrapper()
    RNG = MersenneTwister(1234)
    x = rand(RNG, 1000)
    B = StatUtils.run_bootstrap(1000,sum,x;RNG=RNG)
    @test B.μ ≈ 496.38737230306054
    @test B.σ ≈ 9.052877314367253 
    println("Bootstrapper test passed")

end

function test_fdr()
	pvalues = [0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.00]
	pidx = StatUtils.fdr(pvalues)
	@test pidx == [1,2,3,4]
	println("FDR test passed")
end

function test_groupby()
	groups = StatUtils.groupby([1,2,3], [1,2,3])
	@test Dict{Int64, Array{Int64,1}}(1=>[1], 2=>[2], 3=>[3]) == groups
	println("groupby test passed")
end

@testset "Bootstrapper" begin
    test_bootstrapper()
end
@testset "False Discovery Rate" begin
    test_fdr()
end
@testset "Group by" begin
    test_groupby()
end

@testset "Bootstrap median" begin
    RNG = MersenneTwister(1234)
    x1 = rand(RNG, Float64, (1,100))
    μ,σ = StatUtils.bootstrap_median(x1, 10_000, RNG)
    @test μ[1] ≈ 0.4304271018933901
    @test σ[1] ≈ 0.05446660483997713
end

@testset "Spearman pv" begin
    #no correlation
    RNG = MersenneTwister(1234)
    x1 = rand(RNG, 100)
    x2 = rand(RNG, 100)
    cc,pv = StatUtils.spearmanr(x1,x2, :both, RNG)
    @test cc ≈ -0.0664026402640264
    @test pv ≈ 0.999 

    #small non-significant correlation
    RNG = MersenneTwister(1234)
    x1 = rand(RNG, 100)
    x2 = 0.1*x1 +  0.9*rand(RNG, 100)
    cc,pv = StatUtils.spearmanr(x1,x2, :right, RNG)
    @test cc ≈ 0.031035103510351034
    @test pv ≈ 0.361
    #strong significant correlation
    RNG = MersenneTwister(1234)
    x1 = rand(RNG, 100)
    x2 = 0.5*x1 +  0.5*rand(RNG, 100)
    cc,pv = StatUtils.spearmanr(x1,x2, :right, RNG)
    @test cc ≈ 0.6913171317131713
    @test pv ≈ 0.0
end

@testset "Percentile" begin
    @test_throws ArgumentError Percentile(-0.1)
    @test_throws ArgumentError Percentile(101.0)
    @test typemax(Percentile) == 100.0
    @test typemin(Percentile) == 0.0
end

@testset "Generalized Poisson" begin
    G = GeneralizedPoisson(1.5370664608871476,-0.2658194383776511)
    @test pdf(G,0) ≈ 0.21501092012028916
    @test log(pdf(G,0)) ≈ logpdf(G,0)
    @test log(pdf(G,1)) ≈ logpdf(G,1)
    @test mean(G) ≈ 1.2142857142857142
    @test var(G) ≈ 0.7578397212543557
    @test sum(pdf.(G,0:7)) ≈ 1.0000005196278676
end

@testset "Bootstrap regression" begin
    #generate data
    RNG = MersenneTwister(1234)
    x = range(0.0, stop=1.0, length=20)
    y = 0.5 .+ 0.1*x .+ 0.3*randn(RNG, length(x))
    μ, σ, xx = StatUtils.bootstrap_regression(x,y;RNG=RNG)
    hμ = hash(μ)
    hσ = hash(σ)
    @test hμ == 0x4893fbb4ad08ddbf
    @test hσ == 0x82cc75d75f2c5470
end
