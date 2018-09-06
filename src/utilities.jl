# Assorted utility functions.

"""
    diag(A::Array{T,N}, dim1::Int, dim2::Int)

Collapse `dim1` and `dim2` of `A` into a single dimension containing only the
elements at equal indices.
"""
function diag(A::Array{T,N}, dim1::Int, dim2::Int) where {T,N}
    s = [size(A)...]
    s[dim1] == s[dim2] || throw(DomainError(s, "Dimensions must have the same size."))

    deleteat!(s, dim2)

    result = Array{T}(undef, s...)
    for idx1 in CartesianIndices(A)
        idx1[dim1] == idx1[dim2] || continue
        idx2 = [idx1.I...]
        deleteat!(idx2, dim2)
        result[idx2...] = A[idx1]
    end

    result
end
