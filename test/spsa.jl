using VibronicToolkit: bernoulli, spsa, SpsaInput, SpsaParameters

using Statistics: mean, std

struct SpsaInputTest <: SpsaInput
    x::Float64
    y::Float64
end

Base.:+(si1::SpsaInputTest, si2::SpsaInputTest) = SpsaInputTest(si1.x + si2.x, si1.y + si2.y)
Base.:*(k::Float64, si::SpsaInputTest) = SpsaInputTest(k * si.x, k * si.y)
Base.rand(si::SpsaInputTest) = SpsaInputTest(bernoulli(), bernoulli())

# Rosenbrock function, with a minimum value of 0 at (1, 1).
test_spsa_f(si::SpsaInputTest) = (1 - si.x)^2 + 100 * (si.y - si.x^2)^2

function test_spsa_loss(si::SpsaInputTest)
    # Add noise and sample.
    num_samples = 10_000
    f_samples = [test_spsa_f(si) + randn() for _ in 1:num_samples]
    mean(f_samples), std(f_samples) / sqrt(num_samples)
end

let num_iter = 10_000,
    a = 0.06,
    si = SpsaInputTest(0.5, 3.0),
    (sp, _, _) = SpsaParameters(test_spsa_loss, num_iter, a, si)

    ks = Int[]
    sis = [si]
    opt = spsa(sp, si) do k, si
        push!(ks, k)
        push!(sis, si)
    end

    @test ks == 1:num_iter
    @test length(sis) == num_iter + 1
    @test opt.x == sis[end].x
    @test opt.y == sis[end].y

    # These bounds are generous enough that the tests usually pass, but it's
    # not impossible for them to fail just by chance. If the obtained values
    # are close to the desired ones, these bounds should be increased; if the
    # values are astronomical, the `a` parameter should be decreased.
    @test test_spsa_f(opt) < 5e-2
    @test (opt.x - 1.0)^2 + (opt.y - 1.0)^2 < 5e-1
end
