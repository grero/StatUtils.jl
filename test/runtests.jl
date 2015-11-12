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

test_bootstrapper()
test_fast_bootstrapper()
