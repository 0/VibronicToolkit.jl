# Trial wavefunctions for PIGS.

"""
Trial wavefunction for `S` surfaces and `M` modes.
"""
abstract type TrialWavefunction{S,M} end

"""
Trial wavefunction (for `S` surfaces and `M` modes) that is uniform in space
(i.e. it does not depend on the spatial coordinate).
"""
struct UniformTrialWavefunction{S,M} <: TrialWavefunction{S,M}
    "Surface coefficients of the wavefunction (S)."
    surface_coefs::Vector{Float64}
end

"""
    UniformTrialWavefunction(::System, [surface_coefs::Vector{Float64}])

Create a uniform trial wavefunction with the coefficients `surface_coefs`.

If not provided, `surface_coefs` defaults to all ones.
"""
function UniformTrialWavefunction(::System{S,M}, surface_coefs::Vector{Float64}) where {S,M}
    length(surface_coefs) == S || throw(DomainError(length(surface_coefs), "Incorrect number of coefficients."))
    UniformTrialWavefunction{S,M}(surface_coefs)
end

function UniformTrialWavefunction(sys::System{S,M}) where {S,M}
    UniformTrialWavefunction(sys, ones(S))
end

"""
    trial_spatial(trial::TrialWavefunction, qs::Vector{Float64})

Evaluate `trial` at spatial position `qs`.
"""
function trial_spatial end

function trial_spatial(trial::UniformTrialWavefunction{S,M}, qs::Vector{Float64}) where {S,M}
    length(qs) == M || throw(DomainError(length(qs), "Incorrect number of modes."))
    copy(trial.surface_coefs)
end

"""
    trial_mode(trial::TrialWavefunction{S,M}, basis::Basis{S,M})

Generate a vector for `trial` in `basis`.
"""
function trial_mode end

function trial_mode(trial::UniformTrialWavefunction{S,M}, basis::Basis{S,M}) where {S,M}
    result = zeros(Float64, basis.dim)

    basis_vectors = vectors(basis)
    for idx in 1:basis.dim
        any(isodd, basis_vectors[1:M, idx]) && continue

        k = 1.0
        for m in 1:M
            x = 1.0
            state = basis_vectors[m, idx]
            for j in 1:(state/2)
                x *= (state - j + 1) / (4j)
            end
            k *= sqrt(4pi) * x
        end

        s = basis_vectors[M+1, idx]
        result[idx] = trial.surface_coefs[s] * sqrt(k)
    end

    result
end
