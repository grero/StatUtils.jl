import StatUtils
using Base.Test

function test_bootstrapper()
    srand(1234)
    x = rand(1000)
    B = StatUtils.run_bootstrap(1000,sum,x)
    @test_approx_eq B.μ 497.27583345469264
    @test_approx_eq B.σ 9.145631811573429
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
    x1 = rand(RNG, (1,100))
    μ,σ = StatUtils.bootstrap_median(x1, 10_000, RNG)
    @test μ[1] ≈ 0.42955728184002556
    @test σ[1] ≈ 0.05440210246513024
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
