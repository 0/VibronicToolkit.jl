# Monte Carlo solution for arbitrary systems.

"""
Sampling parameters for a particular system.
"""
abstract type AbstractSamplingParameters{S,M,P} end

"""
Sampling parameters for a particular system at finite temperature.
"""
struct SamplingParameters{S,M,P} <: AbstractSamplingParameters{S,M,P}
    sys::DiagonalSystem{S,M}

    "Imaginary time step."
    tau::Float64

    "Surface weights."
    weights::Weights
    "Precision (inverse covariance) matrices (M, S)."
    precs::Matrix{Matrix{Float64}}
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
    SamplingParameters(sys::DiagonalSystem, beta::Float64, P::Int)

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

    precs = Matrix{Matrix{Float64}}(undef, M, S)
    mvns = Matrix{MvNormal}(undef, M, S)
    for s in 1:S
        for m in 1:M
            mean = sys.ds[m, s] * ones(P)

            # Precision matrix.
            prec = diagm(0 => 2 * Cs[m, s] * ones(P), -1 => -ones(P-1), 1 => -ones(P-1))
            prec[1, end] = -1.0
            prec[end, 1] = -1.0
            prec ./= Ss[m, s]

            # Convariance matrix.
            cov = inv(prec)
            maximum(abs.(cov' - cov)) < 1e-13 || @warn "Asymmetric cov: $(maximum(abs.(cov' - cov)))"
            # Force it to be symmetric.
            cov = (cov + cov') / 2

            precs[m, s] = prec
            mvns[m, s] = MvNormal(mean, cov)
        end
    end

    SamplingParameters{S,M,P}(sys, tau, weights, precs, mvns, Cs, Ss, S_prods)
end

function Base.show(io::IO, sp::SamplingParameters{S,M,P}) where {S,M,P}
    println(io, typeof(sp))
    println(io, "Sampling parameters with $(P) link$(P == 1 ? "" : "s") for a system with $(S) surface$(S == 1 ? "" : "s") and $(M) mode$(M == 1 ? "" : "s").")
    println(io, "Finite temperature.")
    println(io, "imaginary time step (tau): $(sp.tau)")
    println(io, "sampling weights:")
    show_vector(io, weights(sp))
    nothing
end

"""
Sampling parameters for a particular PIGS system.
"""
struct PigsSamplingParameters{S,M,P} <: AbstractSamplingParameters{S,M,P}
    sys::DiagonalSystem{S,M}
    trial::TrialWavefunction{S,M}

    "Imaginary time step."
    tau::Float64

    "Surface weights."
    weights::Weights
    "Precision (inverse covariance) matrices (M, S)."
    precs::Matrix{Matrix{Float64}}
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
    PigsSamplingParameters(sys::DiagonalSystem{S,M}, trial::UniformTrialWavefunction{S,M}, beta::Float64, P::Int)

Generate sampling parameters for `sys` with `trial` propagated by `beta` using
`P` links.
"""
function PigsSamplingParameters(sys::DiagonalSystem{S,M}, trial::UniformTrialWavefunction{S,M}, beta::Float64, P::Int) where {S,M}
    issimple(sys) || throw(DomainError(:sys, "System for sampling must be simple."))

    tau = beta / P

    Zas = abs2.(trial.surface_coefs)
    for s in 1:S
        for m in 1:M
            Zas[s] *= sqrt(2pi / sinh(beta * sys.freq[m, s]))
        end
    end
    preweights = exp.(-beta * (diag(sys.energy) + sys.deltas))
    weights = Weights(preweights .* Zas)

    Cs = cosh.(sys.freq * tau)
    Ss = sinh.(sys.freq * tau)
    S_prods = Float64[prod(Ss[:, s]) for s in 1:S]

    precs = Matrix{Matrix{Float64}}(undef, M, S)
    mvns = Matrix{MvNormal}(undef, M, S)
    for s in 1:S
        for m in 1:M
            mean = sys.ds[m, s] * ones(P+1)

            # Precision matrix.
            prec = diagm(0 => 2 * Cs[m, s] * ones(P+1), -1 => -ones(P), 1 => -ones(P))
            prec[1, 1] /= 2
            prec[end, end] /= 2
            prec ./= Ss[m, s]

            # Convariance matrix.
            cov = inv(prec)
            maximum(abs.(cov' - cov)) < 1e-13 || @warn "Asymmetric cov: $(maximum(abs.(cov' - cov)))"
            # Force it to be symmetric.
            cov = (cov + cov') / 2

            precs[m, s] = prec
            mvns[m, s] = MvNormal(mean, cov)
        end
    end

    PigsSamplingParameters{S,M,P}(sys, trial, tau, weights, precs, mvns, Cs, Ss, S_prods)
end

function Base.show(io::IO, sp::PigsSamplingParameters{S,M,P}) where {S,M,P}
    println(io, typeof(sp))
    println(io, "Sampling parameters with $(P) link$(P == 1 ? "" : "s") for a system with $(S) surface$(S == 1 ? "" : "s") and $(M) mode$(M == 1 ? "" : "s").")
    println(io, "PIGS.")
    println(io, "imaginary time step (tau): $(sp.tau)")
    println(io, "sampling weights:")
    show_vector(io, weights(sp))
    println(io, "trial wavefunction vector: $(sp.trial)")
    nothing
end

"""
    weights(sp::AbstractSamplingParameters)

Compute the normalized surface weights for `sp`.
"""
weights(sp::AbstractSamplingParameters) = sp.weights/sum(sp.weights)

"""
    weights(sys::DiagonalSystem, beta::Float64)

Compute the normalized surface weights for `sys` at `beta` at finite
temperature.
"""
weights(sys::DiagonalSystem, beta::Float64) = weights(SamplingParameters(sys, beta, 2))

"""
Monte Carlo solution for an arbitrary system.
"""
abstract type Sampling <: Solution end

function sampling_matrix_free_particle(sp::AbstractSamplingParameters{S,M,P}, qs1::Vector{Float64}, qs2::Vector{Float64}; scaling=nothing) where {S,M,P}
    FM = zeros(S, S)
    for s in 1:S
        x = -sp.tau * (sp.sys.energy[s, s] + sp.sys.deltas[s])
        for m in 1:M
            q_disp1 = qs1[m] - sp.sys.ds[m, s]
            q_disp2 = qs2[m] - sp.sys.ds[m, s]
            x += -0.5 * ((q_disp1^2 + q_disp2^2) * sp.Cs[m, s] - 2 * q_disp1 * q_disp2) / sp.Ss[m, s]
        end
        FM[s, s] += exp(x) / sqrt(sp.S_prods[s])
    end
    scaling === nothing && (scaling = maximum(FM))
    FM/scaling, scaling
end

function sampling_matrix_interaction(sys::System{S,M}, sp::AbstractSamplingParameters{S,M,P}, qs::Vector{Float64}) where {S,M,P}
    preIM = zeros(S, S)
    for s1 in 1:S
        # Only build the upper triangle.
        for s2 in 1:s1
            if s1 != s2
                preIM[s2, s1] += sys.energy[s2, s1]
                for m in 1:M
                    preIM[s2, s1] += sys.lin[m, s2, s1] * qs[m]
                end
            end
            for m1 in 1:M
                for m2 in 1:M
                    preIM[s2, s1] += 0.5 * sys.quad[m2, m1, s2, s1] * qs[m2] * qs[m1]
                end
            end
        end
    end
    Symmetric(preIM), exp(Symmetric(-sp.tau * preIM))
end

function sampling_matrix_energy(sp::SamplingParameters{S,M,P}, qs1::Vector{Float64}, qs2::Vector{Float64}) where {S,M,P}
    EM2 = zeros(S, S)
    for s in 1:S
        EM2[s, s] -= sp.sys.energy[s, s] + sp.sys.deltas[s]
        for m in 1:M
            q_disp1 = qs1[m] - sp.sys.ds[m, s]
            q_disp2 = qs2[m] - sp.sys.ds[m, s]
            EM2[s, s] -= sp.sys.freq[m, s] * 0.5 * sp.Cs[m, s] / sp.Ss[m, s]
            EM2[s, s] -= sp.sys.freq[m, s] * q_disp1 * q_disp2 * sp.Cs[m, s] / sp.Ss[m, s]^2
            EM2[s, s] += sp.sys.freq[m, s] * 0.5 * (q_disp1^2 + q_disp2^2) / sp.Ss[m, s]^2
        end
    end
    EM2
end

function sampling_matrix_heat_capacity(sp::SamplingParameters{S,M,P}, qs1::Vector{Float64}, qs2::Vector{Float64}) where {S,M,P}
    CM1 = zeros(S, S)
    for s in 1:S
        for m in 1:M
            q_disp1 = qs1[m] - sp.sys.ds[m, s]
            q_disp2 = qs2[m] - sp.sys.ds[m, s]
            CM1[s, s] += sp.sys.freq[m, s]^2 * 0.5 / sp.Ss[m, s]^2
            CM1[s, s] += sp.sys.freq[m, s]^2 * 0.5 * q_disp1 * q_disp2 * (sp.Ss[m, s]^2 + sp.Cs[m, s]^2 + 3) / sp.Ss[m, s]^3
            CM1[s, s] -= sp.sys.freq[m, s]^2 * (q_disp1^2 + q_disp2^2) * sp.Cs[m, s] / sp.Ss[m, s]^3
        end
    end
    CM1
end
