using Distributions
using StatUtils
using Test
using StableRNGs

function test_bootstrapper()
    RNG = StableRNG(1234)
    x = rand(RNG, 1000)
    B = StatUtils.run_bootstrap(1000,sum,x;RNG=RNG)
    @test B.μ ≈ 508.1411134730461
    @test B.σ ≈ 9.045400206380133
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
    RNG = StableRNG(1234)
    x1 = rand(RNG, Float64, (1,100))
    μ,σ = StatUtils.bootstrap_median(x1, 10_000, RNG)
    @test μ[1] ≈ 0.4977072659335834
    @test σ[1] ≈ 0.04304170028498684
end

@testset "Spearman pv" begin
    #no correlation
    RNG = StableRNG(1234)
    x1 = rand(RNG, 100)
    x2 = rand(RNG, 100)
    cc,pv = StatUtils.spearmanr(x1,x2, :both, RNG)
    @test cc ≈ -0.023438343834383438
    @test pv ≈ 0.999

    #small non-significant correlation
    RNG = StableRNG(1234)
    x1 = rand(RNG, 100)
    x2 = 0.1*x1 +  0.9*rand(RNG, 100)
    cc,pv = StatUtils.spearmanr(x1,x2, :right, RNG)
    @test cc ≈ 0.09610561056105611
    @test pv ≈ 0.189
    #strong significant correlation
    RNG = StableRNG(1234)
    x1 = rand(RNG, 100)
    x2 = 0.5*x1 +  0.5*rand(RNG, 100)
    cc,pv = StatUtils.spearmanr(x1,x2, :right, RNG)
    @test cc ≈ 0.6623702370237023
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
    RNG = StableRNG(1234)
    x = range(0.0, stop=1.0, length=50)
    y = 0.5 .+ 0.5*x .+ 0.05*randn(RNG, length(x))
    μ, σ, xx = StatUtils.bootstrap_regression(x,y;RNG=RNG)
    sse = sum(abs2, μ-y)
    sq = sum(abs2, y .- 0.5)
    @test sse/sq < 0.05
end

@testset "Robust regression" begin
    RNG = StableRNG(1234)
    x = range(0.0, stop=1.0, length=20)
    y = 0.5 .+ 0.1*x .+ 0.3*randn(RNG, length(x))
    β0 = rand(RNG, Float64, 2)
    β = StatUtils.robust_regression(x,y, β0)
    @test β[1] ≈ 0.23571370032430217
    @test β[2] ≈ 0.38013900282542895
end


@testset "Angular histogram" begin
    θ = [0.0, π/4, π/2, -π, -π/4, -π/2, -3*π/4]
    counts, binidx, bins = StatUtils.angular_histogram(θ, π/4)
    @test bins ≈ [-π, -3*π/4, -π/2, -π/4, 0.0, π/4, π/2, 3π/4]
    @test counts == [1, 1, 1, 1, 1, 1, 1, 0]
    @test binidx == [5, 6, 7, 1, 4, 3, 2]
end

@testset "Grouping" begin
    x = [1,3,3,1,5,7,7,1,5]
    y = StatUtils.compress(x)
    @test y == [1,2,2,1,3,4,4,1,3]
    X = fill(0.0, 3,4,10,6)
    grouping = fill(0, 10)
    grouping[1:5] .= 1
    grouping[6:10] .= 2
    X[:,:,findall(grouping.==1),:] .= reshape([1.0,2.0,3.0,4.0,5.0], 1,1,5,1)
    X[:,:,findall(grouping.==2),:] .= reshape([6.0, 7.0, 8.0, 9.0, 10.0], 1, 1, 5,1)
    sidx = [7, 10, 2, 5, 8, 6, 1, 3, 9, 4]
    X = X[:,:,sidx,:]
    μ = mean(X, grouping[sidx])
    @test size(μ) == (3,4,2,6)
    @test μ[:,:,1,:] ≈ fill(3.0, 3,4,6)
    @test μ[:,:,2,:] ≈ fill(8.0, 3,4,6)
end
