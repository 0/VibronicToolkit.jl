using VibronicToolkit: jackknife

using Random: MersenneTwister, shuffle
using Statistics: cov, mean, std, var
using StatsBase: Weights

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

# We assume that xs1 and ys1 are not correlated with xs2 and ys2.
function propagate_errors_nonlinear4(f, dfx1, dfx2, dfy1, dfy2, dfx1x1, dfx1y1,
                                     dfx2x2, dfx2y2, dfy1y1, dfy2y2,
                                     xs1, xs2, ys1, ys2)
    mu_x1, mu_x2, mu_y1, mu_y2 = mean(xs1), mean(xs2), mean(ys1), mean(ys2)
    mus = mu_x1, mu_x2, mu_y1, mu_y2
    biased_mean = f(mus...)
    pre_bias1 = dfx1x1(mus...)*var(xs1) + 2dfx1y1(mus...)*cov(xs1,ys1) + dfy1y1(mus...)*var(ys1)
    pre_bias2 = dfx2x2(mus...)*var(xs2) + 2dfx2y2(mus...)*cov(xs2,ys2) + dfy2y2(mus...)*var(ys2)
    unbiased_mean = biased_mean - 0.5*pre_bias1/length(xs1) - 0.5*pre_bias2/length(xs2)
    pre_err1 = dfx1(mus...)^2*var(xs1) + 2dfx1(mus...)*dfy1(mus...)*cov(xs1,ys1) + dfy1(mus...)^2*var(ys1)
    pre_err2 = dfx2(mus...)^2*var(xs2) + 2dfx2(mus...)*dfy2(mus...)*cov(xs2,ys2) + dfy2(mus...)^2*var(ys2)
    err = sqrt(pre_err1/length(xs1) + pre_err2/length(xs2))
    unbiased_mean, err
end

compare_stats((mu1, err1), (mu2, err2)) = isapprox(mu1, mu2; rtol=1e-9) && isapprox(err1, err2; rtol=1e-8)

let rng = MersenneTwister(0x1234),
    xs = randn(rng, 1_000_000),
    ys = -xs .+ 0.5

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

let rng = MersenneTwister(0x1234),
    xs = randn(rng, 1_000_000) .+ 1.0,
    ys = -xs .+ 1.5,
    ws = Weights([3.0]),
    labels = [1 for _ in xs]

    # Nonlinear function of two variables.
    n2a(x, y) = x / y
    dn2ax(x, y) = 1 / y
    dn2ay(x, y) = -x / y^2
    dn2axx(x, y) = 0.0
    dn2axy(x, y) = -1 / y^2
    dn2ayy(x, y) = 2x / y^3
    @test compare_stats(jackknife(n2a, ws, labels, xs, ys),
                        jackknife(n2a, xs, ys))
    @test compare_stats(jackknife(n2a, ws, labels, xs, ys),
                        propagate_errors_nonlinear2(n2a, dn2ax, dn2ay, dn2axx, dn2axy, dn2ayy, xs, ys))
end

let rng = MersenneTwister(0x1234),
    xs1 = randn(rng, 1_000_000),
    ys1 = -xs1 .+ 0.5,
    xs2 = randn(rng, 100_000),
    ys2 = xs2 .+ 1.0,
    ws = Weights([20.0, 1.0]),
    labels = shuffle(rng, [[1 for _ in xs1]; [2 for _ in xs2]])

    xs = Array{Float64}(undef, length(labels))
    ys = Array{Float64}(undef, length(labels))
    idxs = zeros(Int, 2)

    for (i, label) in enumerate(labels)
        idxs[label] += 1
        xs[i] = (label == 1 ? xs1 : xs2)[idxs[label]]
        ys[i] = (label == 1 ? ys1 : ys2)[idxs[label]]
    end

    # Nonlinear function of four variables in linear combinations.
    n4a(x1, x2, y1, y2) = (ws[1] * x1 + ws[2] * x2) / (ws[1] * y1 + ws[2] * y2)
    dn4ax1(x1, x2, y1, y2) = ws[1] / (ws[1] * y1 + ws[2] * y2)
    dn4ax2(x1, x2, y1, y2) = ws[2] / (ws[1] * y1 + ws[2] * y2)
    dn4ay1(x1, x2, y1, y2) = -ws[1] * (ws[1] * x1 + ws[2] * x2) / (ws[1] * y1 + ws[2] * y2)^2
    dn4ay2(x1, x2, y1, y2) = -ws[2] * (ws[1] * x1 + ws[2] * x2) / (ws[1] * y1 + ws[2] * y2)^2
    dn4ax1x1(x1, x2, y1, y2) = 0.0
    dn4ax1y1(x1, x2, y1, y2) = -ws[1]^2 / (ws[1] * y1 + ws[2] * y2)^2
    dn4ax2x2(x1, x2, y1, y2) = 0.0
    dn4ax2y2(x1, x2, y1, y2) = -ws[2]^2 / (ws[1] * y1 + ws[2] * y2)^2
    dn4ay1y1(x1, x2, y1, y2) = 2.0 * ws[1]^2 * (ws[1] * x1 + ws[2] * x2) / (ws[1] * y1 + ws[2] * y2)^3
    dn4ay2y2(x1, x2, y1, y2) = 2.0 * ws[2]^2 * (ws[1] * x1 + ws[2] * x2) / (ws[1] * y1 + ws[2] * y2)^3
    prop_result = propagate_errors_nonlinear4(n4a, dn4ax1, dn4ax2, dn4ay1, dn4ay2,
                                              dn4ax1x1, dn4ax1y1, dn4ax2x2,
                                              dn4ax2y2, dn4ay1y1, dn4ay2y2,
                                              xs1, xs2, ys1, ys2)
    n4a_J1(x, y) = x / y
    @test compare_stats(jackknife(n4a_J1, ws, labels, xs, ys),
                        prop_result)
    n4a_J2(x, ys) = x / ys[1]
    @test compare_stats(jackknife(n4a_J2, ws, labels, xs, [[y] for y in ys]),
                        prop_result)
end

let rng = MersenneTwister(0x1234),
    xs1 = randn(rng, 1_000_000),
    ys1 = -xs1 .+ 0.5,
    # No samples at all for second surface.
    xs2 = Float64[],
    ys2 = Float64[],
    ws = Weights([20.0, 1.0]),
    labels = [1 for _ in xs1]

    # Effective nonlinear function of two variables, since xs2 and ys2 are empty.
    n2a(x, y) = x / y
    dn2ax(x, y) = 1 / y
    dn2ay(x, y) = -x / y^2
    dn2axx(x, y) = 0.0
    dn2axy(x, y) = -1 / y^2
    dn2ayy(x, y) = 2x / y^3
    @test compare_stats(jackknife(n2a, ws, labels, xs1, ys1),
                        propagate_errors_nonlinear2(n2a, dn2ax, dn2ay, dn2axx, dn2axy, dn2ayy, xs1, ys1))
    vn2a(xs) = xs[1] / xs[2]
    @test compare_stats(jackknife(vn2a, ws, labels, [[x, y] for (x, y) in zip(xs1, ys1)]),
                        propagate_errors_nonlinear2(n2a, dn2ax, dn2ay, dn2axx, dn2axy, dn2ayy, xs1, ys1))

    @test_throws DomainError jackknife(n2a, ws, Int[1, 1], xs1, ys1)
    @test_throws DomainError jackknife(n2a, ws, Int[], Float64[], Float64[])
    @test_throws DomainError jackknife(n2a, ws, [3 for _ in xs1], xs1, ys1)
end
