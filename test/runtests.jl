import StatUtils
using Base.Test
using FastAnonymous

function test_bootstrapper()
    srand(1234)
    x = rand(1000)
    B = StatUtils.run_bootstrap(1000,sum,x)
    @test_approx_eq B.μ 497.27583345469264
    @test_approx_eq B.σ 9.145631811573429
    println("Bootstrapper test passed")

end

function test_fast_bootstrapper()
    srand(1234)
    x = rand(1000)

    func = @anon x->begin
       a = 0.0
       for _x in x
         a += _x
       end
       a
       end
    B = StatUtils.run_bootstrap(1000,func,x)

    @test_approx_eq B.μ 497.27583345469264
    @test_approx_eq B.σ 9.145631811573429
    println("Fast bootstrapper test passed")
end

function test_fdr()
	pvalues = [0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.00]
	pidx = StatUtils.fdr(pvalues)
	@test pidx == [1,2,3,4]
	println("FDR test passed")
end

test_bootstrapper()
test_fast_bootstrapper()
test_fdr()
