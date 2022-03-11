function p_norm_sparsity_matrix(M::AbstractArray, ratio::Real, p::Real, min_per_row=0::Integer, min_per_col=0::Integer)
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
    @assert p > 0 "p must be greater than zero"

    n = length(v)

    v = abs.(v)

    nnz_v = sum(v .> 0)

    if nnz_v == 0
        v_pat = spzeros(Int, n)
        return v_pat
    end

    @assert nnz_v >= min_num_nnz "minimum number for non zeros entries bigger than non zero entries"
    @assert nnz_v >= 0 "minimum number for non zeros entries must be greater equal zero"

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
