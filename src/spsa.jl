# Simultaneous perturbation stochastic approximation.

bernoulli() = rand([-1, 1])

"""
Input to a function for SPSA optimization.
"""
abstract type SpsaInput end

Base.:-(si1::SpsaInput, si2::SpsaInput) = si1 + (-1.0 * si2)

"""
Parameters for SPSA optimization.
"""
struct SpsaParameters
    "Number of iterations."
    num_iter::Int

    "Learning rate step offset."
    A::Float64
    "Learning rate coefficient."
    a::Float64
    "Step size coefficient."
    c::Float64

    "Loss function."
    loss::Function
end

"""
    SpsaParameters(loss::Function, num_iter::Int, a::Float64, si::SpsaInput)

Generate appropriate SPSA parameters given the loss function `loss`, number of
iterations `num_iter`, learning rate coefficient `a`, and starting input `si`.

This uses the guidelines in James C. Spall's "implementation" article
(doi:10.1109/7.705889).
"""
function SpsaParameters(loss::Function, num_iter::Int, a::Float64, si::SpsaInput)
    A = div(num_iter, 10)
    l, err = loss(si)

    SpsaParameters(num_iter, A, a, err, loss), l, err
end

"""
    spsa_step(sp::SpsaParameters, k::Int, si::SpsaInput)

Take step `k` of SPSA using `sp` and `si`.

This uses the algorithm in James C. Spall's "implementation" article
(doi:10.1109/7.705889).
"""
function spsa_step(sp::SpsaParameters, k::Int, si::SpsaInput)
    # Gain sequences.
    a_k = sp.a/(sp.A+k)^0.602
    c_k = sp.c/k^0.101

    delta = rand(si)
    l_p, _ = sp.loss(si + c_k * delta)
    l_m, _ = sp.loss(si - c_k * delta)

    si - a_k * (l_p - l_m) / (2c_k) * delta
end

"""
    spsa(f::Function, sp::SpsaParameters, si::SpsaInput)

Optimize `si` using SPSA with `sp`, calling `f` after each iteration.
"""
function spsa(f::Function, sp::SpsaParameters, si::SpsaInput)
    for k in 1:sp.num_iter
        si = spsa_step(sp, k, si)
        f(k, si)
    end
    si
end
