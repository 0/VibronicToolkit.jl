# Solution for a simplified system, ignoring linear coupling between surfaces
# and all quadratic terms.

"""
Exact solution for a simplified system.
"""
immutable Simple <: Solution
    "Partition function."
    Z::Float64
end

"""
    Simple{S,M}(sys::System{S,M}, beta::Beta)

Calculate the solution for the simplified version of `sys` at `beta`.
"""
function Simple{S,M}(sys::System{S,M}, beta::Beta)
    Zs = exp(-beta * diag(sys.energy))

    for s in 1:S
        for m in 1:M
            Zs[s] *= exp(0.5 * beta * sys.lin[m, s, s]^2 / sys.freq[m, s])
            x = exp(-0.5 * beta * sys.freq[m, s])
            Zs[s] *= x / (1-x^2)
        end
    end

    Z = sum(Zs)

    Simple(Z)
end
