using VibronicToolkit: jackknife

using Random: MersenneTwister
using Statistics: cov, mean, std, var

# This uses the formulas from Peter Young's notes
# (https://arxiv.org/abs/1210.3781, https://arxiv.org/pdf/1210.3781v3.pdf).
function propagate_errors_nonlinear1(f, dfx, dfxx, xs)
    mu_x = mean(xs)
    biased_mean = f(mu_x)
    pre_bias = dfxx(mu_x)*var(xs)
    unbiased_mean = biased_mean - 0.5*pre_bias/length(xs)
    pre_err = dfx(mu_x)^2*var(xs)
    err = sqrt(pre_err/length(xs))
    unbiased_mean, err
end

function propagate_errors_nonlinear2(f, dfx, dfy, dfxx, dfxy, dfyy, xs, ys)
    mu_x, mu_y = mean(xs), mean(ys)
    biased_mean = f(mu_x,mu_y)
    pre_bias = dfxx(mu_x,mu_y)*var(xs) + 2dfxy(mu_x,mu_y)*cov(xs,ys) + dfyy(mu_x,mu_y)*var(ys)
    unbiased_mean = biased_mean - 0.5*pre_bias/length(xs)
    pre_err = dfx(mu_x,mu_y)^2*var(xs) + 2dfx(mu_x,mu_y)*dfy(mu_x,mu_y)*cov(xs,ys) + dfy(mu_x,mu_y)^2*var(ys)
    err = sqrt(pre_err/length(xs))
    unbiased_mean, err
end

compare_stats((mu1, err1), (mu2, err2)) = isapprox(mu1, mu2; rtol=1e-11) && isapprox(err1, err2; rtol=1e-6)

let rng = MersenneTwister(0x1234),
    xs = randn(rng, 10000),
    ys = randn(rng, 10000)

    # Linear function of a single variable.
    l1a(x) = 2x-3.7
    @test compare_stats(jackknife(l1a, xs),
                        (l1a(mean(xs)),
                         2std(xs)/sqrt(length(xs))))
    l1b(y) = -1.8y+9
    @test compare_stats(jackknife(l1b, ys),
                        (l1b(mean(ys)),
                         1.8std(ys)/sqrt(length(ys))))

    # Linear function of two variables.
    l2a(x, y) = l1a(x) + l1b(y)
    @test compare_stats(jackknife(l2a, xs, ys),
                        (l2a(mean(xs), mean(ys)),
                         sqrt((2std(xs))^2+(1.8std(ys))^2-2*2*1.8*cov(xs,ys))/sqrt(length(xs))))

    # Nonlinear function of one variable.
    n1a(y) = 1 / l1b(y)
    dn1ay(y) = 1.8 / l1b(y)^2
    dn1ayy(y) = 2*1.8^2 / l1b(y)^3
    @test compare_stats(jackknife(n1a, ys),
                        (propagate_errors_nonlinear1(n1a, dn1ay, dn1ayy, ys)))

    # Nonlinear function of two variables.
    n2a(x, y) = l1a(x) / l1b(y)
    dn2ax(x, y) = 2 / l1b(y)
    dn2ay(x, y) = 1.8*l1a(x) / l1b(y)^2
    dn2axx(x, y) = 0.0
    dn2axy(x, y) = 2*1.8 / l1b(y)^2
    dn2ayy(x, y) = 2*1.8^2*l1a(x) / l1b(y)^3
    @test compare_stats(jackknife(n2a, xs, ys),
                        propagate_errors_nonlinear2(n2a, dn2ax, dn2ay, dn2axx, dn2axy, dn2ayy, xs, ys))

    # As a single vector-valued variable.
    vn2a(xs) = l1a(xs[1]) / l1b(xs[2])
    @test compare_stats(jackknife(vn2a, [[x, y] for (x, y) in zip(xs, ys)]),
                        propagate_errors_nonlinear2(n2a, dn2ax, dn2ay, dn2axx, dn2axy, dn2ayy, xs, ys))

    @test_throws DomainError jackknife(n2a)
    @test_throws DomainError jackknife(n2a, [1.0], [2.0, 3.0])
end
