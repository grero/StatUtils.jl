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


test_bootstrapper()
