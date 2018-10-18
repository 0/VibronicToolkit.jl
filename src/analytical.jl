# Analytical solution for systems with no coupling between surfaces.

"""
Exact solution for an uncoupled system.
"""
abstract type AbstractAnalytical <: Solution end

"""
Exact solution for an uncoupled system at finite temperature.
"""
struct Analytical <: AbstractAnalytical
    "Partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Heat capacity."
    Cv::Float64
end

"""
    Analytical(sys::DiagonalSystem, beta::Float64)

Calculate the solution for `sys` at `beta`.
"""
function Analytical(sys::DiagonalSystem{S,M}, beta::Float64) where {S,M}
    Zs = exp.(-beta * diag(sys.energy))
    E1s = zero(Zs)
    E2s = zero(Zs)

    for s in 1:S
        lambdas, lin_pp = normal_modes(sys, s)
        E_ = sys.energy[s, s]

        for m in 1:M
            x = exp(-0.5 * beta * lambdas[m])
            E_ -= 0.5 * (lin_pp[m] / lambdas[m])^2

            Zs[s] *= exp(0.5 * beta * (lin_pp[m] / lambdas[m])^2)
            Zs[s] *= x / (1-x^2)

            E1s[s] += lambdas[m] * (1+x^2) / (2 * (1-x^2))
        end

        E1s[s] += E_
        E2s[s] += E_^2

        for m1 in 1:M
            x1 = exp(-0.5 * beta * lambdas[m1])

            E2s[s] += 2.0 * E_ * lambdas[m1] * (x1^2+1) / (2 * (1-x1^2))
            E2s[s] += lambdas[m1]^2 * (x1^4+6x1^2+1) / (4 * (1-x1^2)^2)

            for m2 in 1:M
                m1 == m2 && continue

                x2 = exp(-0.5 * beta * lambdas[m2])

                E2s[s] += lambdas[m1] * (x1^2+1) / (2 * (1-x1^2)) * lambdas[m2] * (x2^2+1) / (2 * (1-x2^2))
            end
        end
    end

    Z = sum(Zs)
    E = sum(E1s .* Zs / Z)
    Cv = (sum(E2s .* Zs / Z) - E^2) * beta^2

    Analytical(Z, E, Cv)
end

"""
Exact solution for an uncoupled PIGS system.
"""
struct PigsAnalytical <: AbstractAnalytical
    "Pseudo-partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Von Neumann entanglement entropy."
    SvN::Float64
    "Order-2 Rényi entanglement entropy."
    S2::Float64
end

"""
    PigsAnalytical(sys::DiagonalSystem{S,M}, trial::UniformTrialWavefunction{S,M}, beta::Float64)

Calculate the solution for `sys` with `trial` propagated by `beta`.
"""
function PigsAnalytical(sys::DiagonalSystem{S,M}, trial::UniformTrialWavefunction{S,M}, beta::Float64) where {S,M}
    Zs = exp.(-beta * diag(sys.energy))
    Es = zero(Zs)

    for s in 1:S
        lambdas, lin_pp = normal_modes(sys, s)
        E_ = sys.energy[s, s]

        for m in 1:M
            x = exp(-beta * lambdas[m])
            E_ -= 0.5 * (lin_pp[m] / lambdas[m])^2

            Zs[s] *= exp(0.5 * beta * (lin_pp[m] / lambdas[m])^2)
            Zs[s] *= sqrt(4pi * x / (1-x^2))

            Es[s] += lambdas[m] * (1+x^2) / (2 * (1-x^2))
        end

        Es[s] += E_
    end

    Zs .*= abs2.(trial.surface_coefs)

    rho_surface = trial.surface_coefs * trial.surface_coefs'
    for s1 in 1:S
        for s2 in 1:S
            rho_surface[s2, s1] *= exp(-0.5 * beta * (sys.energy[s1, s1] + sys.energy[s2, s2]))

            lambdas, lin_pp = normal_modes(sys, s1, s2)

            for m in 1:M
                x = exp(-beta * lambdas[m])

                rho_surface[s2, s1] *= exp(0.5 * beta * (lin_pp[m] / lambdas[m])^2)
                rho_surface[s2, s1] *= sqrt(4pi * x / (1-x^2))
            end
        end
    end

    Z = sum(Zs)
    E = sum(Es .* Zs / Z)
    SvN = S_vn(rho_surface / Z)
    S2 = S_renyi(rho_surface / Z)

    PigsAnalytical(Z, E, SvN, S2)
end

"""
    normal_modes(freq::AbstractVector{Float64}, lin::AbstractVector{Float64}, quad::AbstractMatrix{Float64})

Find the normal mode frequencies and linear terms for the coupled oscillators
described by `freq`, `lin`, and `quad`.
"""
function normal_modes(freq::AbstractVector{Float64}, lin::AbstractVector{Float64}, quad::AbstractMatrix{Float64})
    freq_sqrt = sqrt.(freq)
    lin_p = lin .* freq_sqrt
    A = diagm(0 => freq.^2) + quad .* (freq_sqrt * freq_sqrt')
    issymmetric(A) || @warn "Asymmetric A"
    F = eigen(Symmetric(A))
    any(F.values .< 0) && error("Imaginary normal mode frequencies")
    lambdas = sqrt.(F.values)
    T = F.vectors'
    lin_pp = T * lin_p

    lambdas, lin_pp
end

"""
    normal_modes(sys::DiagonalSystem, s::Int)

Find the normal mode frequencies and linear terms for surface `s` in `sys`.
"""
function normal_modes(sys::DiagonalSystem, s::Int)
    normal_modes(sys.freq[:, s], sys.lin[:, s, s], sys.quad[:, :, s, s])
end

"""
    normal_modes(sys::DiagonalSystem, s1::Int, s2::Int)

Find the normal mode frequencies and linear terms for the average of surfaces
`s1` and `s2` in `sys`.
"""
function normal_modes(sys::DiagonalSystem, s1::Int, s2::Int)
    normal_modes(0.5 * (sys.freq[:, s1] + sys.freq[:, s2]),
                 0.5 * (sys.lin[:, s1, s1] + sys.lin[:, s2, s2]),
                 0.5 * (sys.quad[:, :, s1, s1] + sys.quad[:, :, s2, s2]))
end
