# Hamiltonian coefficients.

"""
Arrays of coefficients of various orders for a system with `S` surfaces and `M`
modes.
"""
struct HamiltonianCoefficients{S,M}
    """
    Dictionary of coefficients. Each key is a non-negative integer. Each value
    has the shape (M, ..., M, S, S) with the number of mode dimensions equal to
    the key.
    """
    d::Dict{Int,Array{Float64}}
    "List of orders that are populated entirely with zeros."
    zerolist::Vector{Int}

    function HamiltonianCoefficients{S,M}(d::Dict{Int,Array{Float64}}) where {S,M}
        check_shape(S, M, d)
        zerolist = [ord for (ord, val) in d if iszero(val)]
        new{S,M}(d, zerolist)
    end

    function HamiltonianCoefficients(d::Dict{Int,Array{Float64}})
        S, M = get_shape(d)
        HamiltonianCoefficients{S,M}(d)
    end
end

function HamiltonianCoefficients{S,M}(As::AbstractArray{Float64}...) where {S,M}
    HamiltonianCoefficients{S,M}(Dict{Int,Array{Float64}}(zip(0:(length(As)-1), As)))
end

function HamiltonianCoefficients(As::AbstractArray{Float64}...)
    HamiltonianCoefficients(Dict{Int,Array{Float64}}(zip(0:(length(As)-1), As)))
end

function Base.haskey(coef::HamiltonianCoefficients, ord::Int)
    haskey(coef.d, ord) || return false
    ord in coef.zerolist && return false
    true
end

function Base.getindex(coef::HamiltonianCoefficients{S,M}, ord::Int) where {S,M}
    haskey(coef.d, ord) && return coef.d[ord]
    # If we want the array for an order that doesn't exist, default to zeros.
    zeros(Float64, [M for _ in 1:ord]..., S, S)
end

Base.keys(coef::HamiltonianCoefficients) = sort([ord for ord in keys(coef.d) if haskey(coef, ord)])

function Base.iterate(coef::HamiltonianCoefficients, state=1)
    ords = keys(coef)
    state > length(ords) && return nothing
    ord = ords[state]
    ord => coef[ord], state+1
end

Base.eltype(::Type{HamiltonianCoefficients}) = Array{Float64}
Base.length(coef::HamiltonianCoefficients) = length(keys(coef))

function Base.:*(x::Float64, coef::HamiltonianCoefficients{S,M}) where {S,M}
    HamiltonianCoefficients{S,M}(Dict(ord => x * val for (ord, val) in coef.d))
end

"""
    get_shape(d::Dict{Int,Array{Float64}})

Extract the system shape from the coefficient arrays in `d`.
"""
function get_shape(d::Dict{Int,Array{Float64}})
    ord_max = maximum(keys(d))
    ord_max >= 1 || throw(DomainError(keys(d), "Must supply at least first order coefficients."))
    M, S = size(d[ord_max])[1], size(d[ord_max])[end]
    S, M
end

"""
    check_shape(S::Int, M::Int, d::Dict{Int,Array{Float64}})

Throw the appropriate exception if the coefficient arrays in `d` do not have
`S` surfaces and `M` modes.
"""
function check_shape(S::Int, M::Int, d::Dict{Int,Array{Float64}})
    S >= 1 || throw(DomainError(S, "At least 1 surface."))
    M >= 1 || throw(DomainError(M, "At least 1 mode."))
    for (ord, val) in d
        size(val) == tuple([M for _ in 1:ord]..., S, S) || throw(DomainError(size(val), "Unexpected coefficient dimensions."))
    end
    nothing
end

"""
    isdiag(coef::HamiltonianCoefficients)

Determine whether the coefficients in `coef` are diagonal in surfaces.
"""
function isdiag(coef::HamiltonianCoefficients)
    for (ord, val) in coef
        for idx in mode_indices(val)
            isdiag(val[idx, :, :]) || return false
        end
    end
    true
end

"""
    mode_indices(xs::Array{Float64})

Find the `CartesianIndices` corresponding to the mode dimensions of `xs` (i.e.
all but the last 2 dimensions).
"""
mode_indices(xs::Array{Float64}) = CartesianIndices(size(xs)[1:end-2])

"""
    coefficients(f::Function, coef::HamiltonianCoefficients)

Call `f` for every non-zero coefficient in `coef`.
"""
function coefficients(f::Function, coef::HamiltonianCoefficients{S,M}) where {S,M}
    for (ord, val) in coef
        for s1 in 1:S
            for s2 in 1:S
                for idx in mode_indices(val)
                    val[idx, s2, s1] == 0 && continue
                    f(ord, s1, s2, Tuple(idx), val[idx, s2, s1])
                end
            end
        end
    end
end
