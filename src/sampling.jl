# Monte Carlo solution for arbitrary systems.

"""
Sampling parameters for a particular system.
"""
struct SamplingParameters{S,M,P}
    sys::DiagonalSystem{S,M}

    "Imaginary time step."
    tau::Float64

    "Surface weights."
    weights::Weights
    "Multivariate normal distributions (M, S)."
    mvns::Matrix{MvNormal}

    "Precomputed cosh factors (M, S)."
    Cs::Matrix{Float64}
    "Precomputed sinh factors (M, S)."
    Ss::Matrix{Float64}
    "Precomputed products of sinh factors (S)."
    S_prods::Vector{Float64}
end

"""
    SamplingParameters(sys::DiagonalSystem{S,M}, beta::Float64, P::Int)

Generate sampling parameters for `sys` at `beta` with `P` links.
"""
function SamplingParameters(sys::DiagonalSystem{S,M}, beta::Float64, P::Int) where {S,M}
    issimple(sys) || throw(DomainError(:sys, "System for sampling must be simple."))

    tau = beta / P

    Zas = ones(S)
    for s in 1:S
        for m in 1:M
            Zas[s] *= 1.0 / (2*sinh(0.5 * beta * sys.freq[m, s]))
        end
    end
    preweights = exp.(-beta * (diag(sys.energy) + sys.deltas))
    weights = Weights(preweights .* Zas)

    Cs = cosh.(sys.freq * tau)
    Ss = sinh.(sys.freq * tau)
    S_prods = Float64[prod(Ss[:, s]) for s in 1:S]

    mvns = Matrix{MvNormal}(undef, M, S)
    for s in 1:S
        for m in 1:M
            mean = sys.ds[m, s] * ones(P)

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

    SamplingParameters{S,M,P}(sys, tau, weights, mvns, Cs, Ss, S_prods)
end

function Base.show(io::IO, sp::SamplingParameters{S,M,P}) where {S,M,P}
    println(io, typeof(sp))
    println(io, "Sampling parameters with $(P) bead$(P == 1 ? "" : "s") for a system with $(S) surface$(S == 1 ? "" : "s") and $(M) mode$(M == 1 ? "" : "s").")
    println(io, "imaginary time step (tau): $(sp.tau)")
    println(io, "sampling weights:")
    show_vector(io, sp.weights/sum(sp.weights))
    nothing
end

"""
Monte Carlo solution for an arbitrary system.
"""
abstract type Sampling <: Solution end
