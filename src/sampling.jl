# Monte Carlo solution for arbitrary systems.

"""
Sampling parameters for a particular system.
"""
struct SamplingParameters{S,M,P}
    sys::System{S,M}

    "Imaginary time step."
    tau::Float64

    "Surface weights."
    weights::Weights
    "Multivariate normal distributions."
    mvns::Matrix{MvNormal}

    "Energy offsets due to linear terms."
    deltas::Vector{Float64}
    "Position offsets due to linear terms."
    ds::Matrix{Float64}

    "Precomputed cosh factors."
    Cs::Matrix{Float64}
    "Precomputed sinh factors."
    Ss::Matrix{Float64}
    "Precomputed products of sinh factors."
    S_prods::Vector{Float64}
end

"""
    SamplingParameters(sys::System{S,M}, beta::Float64, P::Int)

Generate sampling parameters for `sys` at `beta` with `P` links.
"""
function SamplingParameters(sys::System{S,M}, beta::Float64, P::Int) where {S,M}
    tau = beta / P

    Zas = ones(S)
    deltas = zeros(S)
    ds = zeros(M, S)
    for s in 1:S
        for m in 1:M
            Zas[s] *= 1.0 / (2*sinh(0.5 * beta * sys.freq[m, s]))
            deltas[s] += -0.5 * sys.lin[m, s, s].^2 / sys.freq[m, s]
            ds[m, s] = -sys.lin[m, s, s] / sys.freq[m, s]
        end
    end

    preweights = exp.(-beta * (diag(sys.energy) + deltas))
    weights = Weights(preweights .* Zas)

    Cs = cosh.(sys.freq * tau)
    Ss = sinh.(sys.freq * tau)
    S_prods = Float64[prod(Ss[:, s]) for s in 1:S]

    mvns = Matrix{MvNormal}(undef, M, S)
    for s in 1:S
        for m in 1:M
            mean = ds[m, s] * ones(P)

            # Precision matrix.
            prec = diagm(0 => 2 * Cs[m, s] * ones(P), -1 => -ones(P-1), 1 => -ones(P-1))
            prec[1, end] = -1.0
            prec[end, 1] = -1.0

            # Convariance matrix.
            cov = inv(prec / Ss[m, s])
            maximum(abs.(cov' - cov)) < 1e-13 || @warn "Asymmetric cov: $(maximum(abs.(cov' - cov)))"
            # Force it to be symmetric.
            cov = (cov + cov') / 2

            mvns[m, s] = MvNormal(mean, cov)
        end
    end

    SamplingParameters{S,M,P}(sys, tau, weights, mvns, deltas, ds, Cs, Ss, S_prods)
end

"""
Monte Carlo solution for an arbitrary system.
"""
abstract type Sampling <: Solution end
