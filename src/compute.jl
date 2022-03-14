"""
    ssa_compute(M::AbstractArray, ratio::Real, p::Real, max_num_bins::Integer, 
        impose_null_spaces=false::Bool)

Compute a SparseMatrixCSC{Float64, Int64} `X` for a matrix `M`, that is sparse and spectrally close to M,
especially in the lower end of the singular value spectrum.

Wrapper function over all the sparse spectral approximation (ssa) functionality.
Compute a sparse approximation for a matrix `M`. First the sparsity pattern is computetd with a
 `ratio` in [0,1] that determines the sparsity (0 means more sparse, 1 means less) 
and with `p` in (0,Inf] that determines which norm is used to find the sparsity pattern.
Next the sparsity pattern is used to find a binning pattern, the non-negative integer `max_num_bins`
determines whats the maximum number of bins (200-1000 is a reasonable choice). `max_num_bins == 0`
means no binning is performed and can lead to significant slowdown.
Then an optimization problem, with constraints given by the binning pattern, is solved to find an 
sparse approximation `X` for `M`.
If `impose_null_spaces == true`, another optimization problem is solved that ensures that `X` and
`M` have the same null space.

See also: [`p_norm_sparsity_matrix`](@ref), [`bin_sparse_matrix!`](@ref).

# Examples
```jldoctest
julia> ssa_compute([16.99 65; 0.1 17.01], 0.6, 2, 200)
2×2 SparseMatrixCSC{Float64, Int64} with 3 stored entries:
 16.8041  64.2499
   ⋅      16.8041
```
"""
function ssa_compute(M::AbstractArray{T}, ratio::Real, p::Real, max_num_bins::Integer, 
    impose_null_spaces=false::Bool) where{T}

    if T <: Integer
        M = 1.0*M
    end
    
    pinv_M, rnull, lnull = pinv_qr(M)

    # sparsity pattern
    num_near_zero_rows, num_near_zero_cols = near_zero_row_col(M)
    min_per_row = max(0, min(size(rnull)[2] - num_near_zero_cols, size(M, 2)))
    min_per_col = max(0, min(size(lnull)[2] - num_near_zero_rows, size(M, 1)))
    M_id = p_norm_sparsity_matrix(M, ratio, p, min_per_row, min_per_col)

    # binning pattern
    bin_sparse_matrix!(M, M_id, max_num_bins)

    # minimization
    X = ssa_minimization(M, M_id, pinv_M)
    impose_null_spaces && ssa_impose_action!(X, M, rnull, lnull)
    return X
end


function near_zero_row_col(M::AbstractArray{T}) where {T}
    max_in_col = maximum(M, dims=1)
    max_in_row = maximum(M, dims=2)
    A_max = maximum(max_in_col)

    tol = A_max * eps(T) * 100  # magic constant from matlab

    num_near_zero_cols = sum(max_in_col .<= tol)
    num_near_zero_rows = sum(max_in_row .<= tol)
    return num_near_zero_rows, num_near_zero_cols
end


function ssa_minimization!(M, M_id, pinv_M)
    pinv_MTM = pinv_M * pinv_M'
    pinv_MMT = pinv_M' * pinv_M
    LS_M, LS_b = ssa_system_no_null(M, M_id, pinv_MTM, pinv_MMT)
    LS_M, = pinv_qr(LS_M)
    LS_x = LS_M * LS_b

    M[:, :] = ssa_unknown_to_matrix(LS_x, M_id)
    return M
end


function ssa_minimization(M, M_id, pinv_M)
    pinv_MTM = pinv_M * pinv_M'
    pinv_MMT = pinv_M' * pinv_M
    LS_M, LS_b = ssa_system_no_null(M, M_id, pinv_MTM, pinv_MMT)
    LS_M, = pinv_qr(LS_M) 

    LS_x = LS_M * LS_b

    Y = ssa_unknown_to_matrix(LS_x, M_id)
    return Y
end


function ssa_system_no_null(M::AbstractArray, M_id::AbstractSparseMatrixCSC, pinv_MTM::AbstractMatrix, pinv_MMT::AbstractMatrix)

    N = length(unique(M_id.nzval))

    D = M * pinv_MTM + pinv_MMT * M

    LS_A = zeros(Float64, (N, N))
    LS_b = zeros(Float64, N)

    for idx in findall(!iszero, M_id)
        J, L = findnz(M_id[Tuple(idx)[1], :])
        k = 1
        for l in L
            @views LS_A[M_id[idx], l] = LS_A[M_id[idx], l] + pinv_MTM[J[k], Tuple(idx)[2]]
            k += 1
        end

        J, L = findnz(M_id[:, Tuple(idx)[2]])
        k = 1
        for l in L
            @views LS_A[M_id[idx], l] = LS_A[M_id[idx], l] + pinv_MMT[Tuple(idx)[1], J[k]]
            k += 1
        end
    end

    @views for i in 1:N
        LS_b[i] = sum(D[M_id .== i])
    end

    return LS_A, LS_b
end


function ssa_unknown_to_matrix(x::AbstractVector, M_id::AbstractSparseMatrixCSC)
    m, n = size(M_id)
    M_id_i, M_id_j, M_id_v = findnz(M_id)
    return sparse(M_id_i, M_id_j, x[M_id_v], m, n)
end


function ssa_impose_action!(Y::AbstractSparseMatrixCSC, M::AbstractArray, R_mat::AbstractMatrix, L_mat::AbstractMatrix,
     tol=eps(Float64)::Real, max_iters=1000::Integer)
    m, n = size(Y)
    _, nR = size(R_mat)
    _, nL = size(L_mat)

    # Quantities don't change when iterating
    norm_R_mat = norm(R_mat)
    norm_L_mat = norm(L_mat)

    M_R_mat = M * R_mat
    M_L_mat = M' * L_mat

    nnz_Y = count(!iszero, Y)
    Y_i, Y_j, _ = findnz(Y)

    # Quantities that do change when iterating
    Lag_R = zeros(m, nR)
    Lag_L = zeros(n, nL)

    d_R = Y * R_mat - M_R_mat
    d_L = Y' * L_mat - M_L_mat

    q_R = - d_R
    q_L = - d_L

    d_R_proj_vec = zeros(nnz_Y, 1)
    d_L_proj_vec = zeros(nnz_Y, 1)
    d_proj_vec = zeros(nnz_Y, 1)

    iters = 1
    for i in 1:max_iters
        if norm(d_R) <= tol * n * norm(Y - M) * norm_R_mat  norm(d_L) <= tol * m * norm(Y - M) * norm_L_mat
            break
        end

        @views for k in 1:nnz_Y
            d_R_proj_vec[k] = R_mat[Y_j[k], :]' * d_R[Y_i[k], :]
            d_L_proj_vec[k] = L_mat[Y_i[k], :]' * d_L[Y_j[k], :]
        end

        d_proj_vec = d_R_proj_vec .+ d_L_proj_vec

        alpha = (norm(q_R)^2 + norm(q_L)^2) / norm(d_proj_vec)^2

        Lag_R = Lag_R .+ alpha .* d_R
        Lag_L = Lag_L .+ alpha .* d_L

        @views for k in 1:nnz_Y
            Y[Y_i[k], Y_j[k]] = Y[Y_i[k], Y_j[k]] - alpha * d_proj_vec[k]
        end

        norm2_q_RL_old = norm(q_R)^2 + norm(q_L)^2

        q_R = M_R_mat .- Y * R_mat
        q_L = M_L_mat .- Y' * L_mat

        norm2_q_RL_new = norm(q_R)^2 + norm(q_L)^2

        beta = norm2_q_RL_new / norm2_q_RL_old
        d_R = beta .* d_R .- q_R
        d_L = beta .* d_L .- q_L

        iters = iters + 1
    end
    return Y
end
