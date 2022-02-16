function p_norm_sparsity_matrix(M::AbstractMatrix, ratio::Real, p::Real, min_per_row::Integer, min_per_col::Integer)
    # todo implement all different cases
    # todo make almost symmetric matrices symmetric
    m, n = size(M)
    M_pat = spzeros(Int, m, n)
    for i in 1:m
        M_pat[i, :] = p_norm_sparsity_vector(M[i, :], ratio, p, min_per_row)
    end
    for j in 1:n
        M_pat[:, j] = broadcast(max, p_norm_sparsity_vector(M[:, j], ratio, p, min_per_col), M_pat[:, j])
    end
    return M_pat
end


function p_norm_sparsity_vector(v::AbstractVector, ratio::Real, p::Real, min_num_nnz::Integer)
    # todo implement all different cases
    # todo p = inf, p < 1
    n = length(v)

    v = broadcast(abs, v)
    v = v.^p

    nnz_v = sum(v .> 0) # todo if nnz_v == 0

    idx = sortperm(v)
    v = v[idx]

    csum = cumsum(v)
    sum_threshold = csum[end] * (1 - ratio) ^ max(1, p)

    ub = upper_bound(csum, sum_threshold)
    num_needed = max(min_num_nnz, n - ub + 1)
    data = ones(Int, num_needed)
    v_pat = sparsevec(idx[n - num_needed + 1:n], data, n)
    return v_pat
end


function upper_bound(v::AbstractVector, value::Real)
    n = length(v)
    low = 0
    high = n-1
    while low <= high
        mid = floor(Int, (low + high) / 2)
        if value < v[mid]
            high = mid - 1
        else
            low = mid + 1
        end
    end
    return low
end
