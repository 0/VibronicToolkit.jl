# Utilities for entanglement entropy.

"""
    S_vn(eigvals::AbstractVector{Float64})

Von Neumann entropy of `eigvals`.
"""
S_vn(eigvals::AbstractVector{Float64}) = -sum(x * log(x) for x in eigvals if x > 0)

"""
    S_vn(rho::AbstractMatrix{Float64})

Von Neumann entropy of `rho`.
"""
S_vn(rho::AbstractMatrix{Float64}) = rho |> eigvals |> S_vn

"""
    S_renyi(eigvals::AbstractVector{Float64}, alpha=2)

Order-`alpha` Rényi entropy of `eigvals`.
"""
S_renyi(eigvals::AbstractVector{Float64}; alpha=2) = log(sum(eigvals.^alpha)) / (1 - alpha)

"""
    S_renyi(rho::AbstractMatrix{Float64}; alpha=2)

Order-`alpha` Rényi entropy of `rho`.
"""
S_renyi(rho::AbstractMatrix{Float64}; alpha=2) = S_renyi(rho |> eigvals; alpha=alpha)
