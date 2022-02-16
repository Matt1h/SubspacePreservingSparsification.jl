function ssa_compute(M::AbstractMatrix, ratio::Real, p::Real, max_num_bins::Integer, impose_null_spaces=false::Bool)
    pinv_M, rnull, lnull = pinv_qr(M)

    # sparsity pattern
    num_near_zero_rows, num_near_zero_cols = near_zero_row_col(M)
    min_per_row = max(0, min(size(rnull)[2] - num_near_zero_cols, size(M)[2]))
    min_per_col = max(0, min(size(lnull)[2] - num_near_zero_rows, size(M)[1]))
    M_id = p_norm_sparsity_matrix(M, ratio, p, min_per_row, min_per_col)

    # binning pattern
    bin_sparse_matrix!(M, M_id, max_num_bins) 

    # minimization
    if impose_null_spaces
        X = ssa_minimization(M, M_id, pinv_M)
        ssa_impose_action!(X, M, rnull, lnull)
    else
        X = ssa_minimization(M, M_id, pinv_M)
        return X
    end
    return X
end


function near_zero_row_col(M::AbstractMatrix{T}) where {T}
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


function ssa_system_no_null(M::AbstractMatrix, M_id::AbstractSparseMatrixCSC, pinv_MTM::AbstractMatrix, pinv_MMT::AbstractMatrix)
    M_id_i, M_id_j, M_id_v = findnz(M_id)  # TODO: compute this only one time?

    M_id_uniq = unique(M_id_v)
    N = length(M_id_uniq)

    # number of most frequently occurring dof id
    max_id_occurrences = sum(M_id_v .== mode(M_id_v))

    M_id_loc_row = zeros(Int, N, max_id_occurrences)
    M_id_loc_col = zeros(Int, N, max_id_occurrences)
    M_id_loc_next = ones(Int, N)

    M_id_nnz = count(!iszero, M_id)
    
    for k in 1:M_id_nnz
        M_id_val = M_id_v[k]
        M_id_loc_row[M_id_val, M_id_loc_next[M_id_val]] = M_id_i[k]
        M_id_loc_col[M_id_val, M_id_loc_next[M_id_val]] = M_id_j[k]
        M_id_loc_next[M_id_val] = M_id_loc_next[M_id_val] + 1
    end
    
    M_id_count = M_id_loc_next .- 1

    M__pinv_ATA = M * pinv_MTM
    pinv_MMT__A = pinv_MMT * M

    LS_A = Matrix{Float64}(undef, N, N)
    LS_b = Vector{Float64}(undef, N)

    for i in 1:N
        for j in 1:N

            tmp = 0
            for K in 1:M_id_count[i]
                for L in 1:M_id_count[j]

                    u_row = M_id_loc_row[i, K]
                    u_col = M_id_loc_col[i, K]

                    v_row = M_id_loc_row[j, L]
                    v_col = M_id_loc_col[j, L]

                    if u_col == v_col
                        tmp = tmp + pinv_MMT[u_row, v_row]
                    end
                    if u_row == v_row
                        tmp = tmp + pinv_MTM[u_col, v_col]
                    end
                end
            end
            LS_A[i, j] = tmp
        end
    end


    for i in 1:N
        tmp = 0
        for K in 1:M_id_count[i]
            row = M_id_loc_row[i, K]
            col = M_id_loc_col[i, K]

            tmp = tmp + M__pinv_ATA[row, col] + pinv_MMT__A[row, col]
        end
        LS_b[i] = tmp
    end

    return LS_A, LS_b
end


function ssa_unknown_to_matrix(x::AbstractVector, M_id::AbstractSparseMatrixCSC)
    M_id_i, M_id_j, M_id_v = findnz(M_id)
    return sparse(M_id_i, M_id_j, x[M_id_v])
end


function ssa_impose_action!(Y::AbstractSparseMatrixCSC, M::AbstractMatrix, R_mat::AbstractMatrix, L_mat::AbstractMatrix,
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
