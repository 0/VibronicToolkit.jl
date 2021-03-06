# Analytical solution for systems with no coupling between surfaces.

"""
Exact solution for an uncoupled system.
"""
struct Analytical <: Solution
    "Partition function."
    Z::Float64
    "Energy."
    E::Float64
    "Heat capacity."
    Cv::Float64
end

"""
    Analytical(sys::DiagonalSystem{S,M}, beta::Float64)

Calculate the solution for `sys` at `beta`.
"""
function Analytical(sys::DiagonalSystem{S,M}, beta::Float64) where {S,M}
    Zs = exp.(-beta * diag(sys.energy))
    E1s = zero(Zs)
    E2s = zero(Zs)

    for s in 1:S
        # Find the normal mode frequencies for this surface.
        freq_sqrt = sqrt.(sys.freq[:, s])
        lin_p = sys.lin[:, s, s] .* freq_sqrt
        A = diagm(0 => sys.freq[:, s].^2) + sys.quad[:, :, s, s] .* (freq_sqrt * freq_sqrt')
        issymmetric(A) || @warn "Asymmetric A"
        F = eigen(Symmetric(A))
        any(F.values .< 0) && error("Imaginary normal mode frequencies")
        lambdas = sqrt.(F.values)
        T = F.vectors'
        lin_pp = T * lin_p

        E_ = sys.energy[s, s]

        for m in 1:M
            x = exp(-0.5 * beta * lambdas[m])
            E_ -= 0.5 * (lin_pp[m] / lambdas[m])^2

            Zs[s] *= exp(0.5 * beta * (lin_pp[m] / lambdas[m])^2)
            Zs[s] *= x / (1-x^2)

            E1s[s] += lambdas[m] * (x^2+1) / (2 * (1-x^2))
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
