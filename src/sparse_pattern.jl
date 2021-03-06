"""
    sparsity_pattern(M::AbstractArray, ratio::Real, p::Real, min_per_row=0::Integer, min_per_col=0::Integer)

Compute a pattern matrix `M_pat` which is a p-norm sparsity pattern for `M`.

`M_pat` is a sparse `SparseMatrixCSC{Int64, Int64}` of the same shape as `M` and
contains only 0 or 1.

`ratio` is a measure for how sparse the pattern is supposed to be and must be in ``[0,1]``.  
`p` determines which norm is used and should be in ``(0, \\infty ]``. `min_per_row` is the minimum number 
of non-zeros needed per row and `min_per_col` is the minimum number of non-zeros needed per column, 
Defaults are `0`.

See also: [`sparsify`](@ref), [`bin_sparse_matrix!`](@ref).

# Examples
```julia-repl
julia> sparsity_pattern([4 1; 0.1 17], 0.6, 2)
2×2 SparseMatrixCSC{Int64, Int64} with 2 stored entries:
 1  ⋅
 ⋅  1
```
"""
function sparsity_pattern(M::AbstractArray, ratio::Real, p::Real, min_per_row=0::Integer, min_per_col=0::Integer)
    0 <= ratio <= 1 || throw(DomainError(ratio, "ratio must be in [0,1]"))
    p > 0 || throw(DomainError(p, "p must be greater than zero"))
    min_per_row >= 0|| throw(DomainError("min_per_row must be greater equal zero"))
    min_per_col >= 0|| throw(DomainError("min_per_col must be greater equal zero"))    
    
    # todo make almost symmetric matrices symmetric
    m = size(M, 1)
    n = size(M, 2)

    M_pat = spzeros(Int, m, n)
    @views for i in 1:m
        M_pat[i, :] = p_norm_sparsity_vector(M[i, :], ratio, p, min_per_row)
    end
    @views for j in 1:n
        M_pat[:, j] = max.(p_norm_sparsity_vector(M[:, j], ratio, p, min_per_col), M_pat[:, j])
    end
    return M_pat
end


function p_norm_sparsity_vector(v::AbstractVector, ratio::Real, p::Real, min_num_nnz::Integer)
    0 <= ratio <= 1 || throw(DomainError(ratio, "ratio must be in [0,1]"))
    p > 0 || throw(DomainError(p, "p must be greater than zero"))
    min_num_nnz >= 0|| throw(DomainError("min_num_nnz must be greater equal zero"))

    n = length(v)

    v = abs.(v)

    nnz_v = sum(v .> 0)

    if nnz_v == 0
        v_pat = spzeros(Int, n)
        return v_pat
    end

    nnz_v >= min_num_nnz || throw(ArgumentError("minimum number for non zeros entries bigger than non zero entries in v"))

    idx = sortperm(v)
    v = v[idx]

    if p == Inf
        if min_num_nnz == nnz_v 
            num_needed = min_num_nnz
        else
            discarded_bound = (1-ratio)*v[end]
            num_needed = n - upper_bound(v[1:end-min_num_nnz], discarded_bound) + 1
        end
    else
        v = v.^p
        csum = cumsum(v)
        sum_threshold = csum[end] * (1 - ratio) ^ max(1, p)

        ub = upper_bound(csum, sum_threshold)
        num_needed = max(min_num_nnz, n - ub + 1)
    end

    data = ones(Int, num_needed)
    v_pat = sparsevec(idx[n - num_needed + 1:n], data, n)
    return v_pat
end


function upper_bound(v::AbstractVector, value::Real)
    n = length(v)

    if n == 0
        low = 0
    else
        low = 1
        high = n
        while low <= high
            mid = floor(Int, (low + high) / 2)
            if value < v[mid]
                high = mid - 1
            else
                low = mid + 1
            end
        end
    end
    return low
end
