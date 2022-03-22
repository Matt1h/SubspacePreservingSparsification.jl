"""
    sparsify(M::AbstractArray, ratio::Real, p::Real, max_num_bins::Integer,
        impose_null_spaces=false::Bool)

Compute a `SparseMatrixCSC{Float64, Int64}` `X` for a matrix `M`, that is sparse and spectrally close to `M`,
especially in the lower end of the singular value spectrum.

Wrapper function over all functionality.

Compute a sparse approximation for a matrix `M`. First the sparsity pattern is computed with a
`ratio` in ``[0,1]`` that determines the sparsity (`0` means more sparse, `1` means less)
and with `p` in ``(0, \\infty ]`` that determines which norm is used to find the sparsity pattern.
Next the sparsity pattern is used to find a binning pattern, the non-negative integer `max_num_bins`
determines whats the maximum number of bins is (200-1000 is a reasonable choice). `max_num_bins=0`
means no binning is performed and can lead to significant slowdown.
Then an optimization problem, with constraints given by the binning pattern, is solved to find a
sparse approximation `X` for `M`.
If `impose_null_spaces=true`, another optimization problem is solved that ensures that `X` and
`M` have the same null space.

See also: [`sparsity_pattern`](@ref), [`bin_sparse_matrix!`](@ref).

# Examples
```julia-repl
julia> sparsify([16.99 65; 0.1 17.01], 0.6, 2, 200)
2×2 SparseMatrixCSC{Float64, Int64} with 3 stored entries:
 16.8041  64.2499
   ⋅      16.8041
```
"""
function sparsify(M::AbstractArray{T}, ratio::Real, p::Real, max_num_bins::Integer, 
    impose_null_spaces=false::Bool) where{T}
    
    pinv_M, rnull, lnull = pinv_qr(M)

    # sparsity pattern
    num_near_zero_rows, num_near_zero_cols = near_zero_row_col(M)
    min_per_row = max(0, min(size(rnull)[2] - num_near_zero_cols, size(M, 2)))
    min_per_col = max(0, min(size(lnull)[2] - num_near_zero_rows, size(M, 1)))

    M_id = sparsity_pattern(M, ratio, p, min_per_row, min_per_col)

    # binning pattern
    bin_sparse_matrix!(M, M_id, max_num_bins)

    # minimization
    pinv_MTM = pinv_M * pinv_M'
    pinv_MMT = pinv_M' * pinv_M
    D = M * pinv_MTM + pinv_MMT * M
    X = binned_minimization(M, M_id, pinv_MTM, pinv_MMT, D)
    impose_null_spaces && sps_impose_action!(X, M, rnull, lnull)
    return X
end


function sparsify(M::AbstractArray{<:Int}, ratio::Real, p::Real, max_num_bins::Integer, impose_null_spaces=false::Bool)
    return sparsify(float(M), ratio, p, max_num_bins, impose_null_spaces)
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


function near_zero_row_col(M::AbstractArray{<:Int})
    near_zero_row_col(float(M))
end


function binned_minimization(M::AbstractArray, M_id::AbstractSparseMatrixCSC,
    B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix)
    m = size(M, 1)
    n = size(M, 2)

    size(M_id, 1) == size(M, 1) && size(M_id, 2) == size(M, 2) || 
    throw(ArgumentError("M_id must be the same size as M"))
    0 <= minimum(M_id) || throw(DomainError(minimum(M_id), "M_id contains negative integers"))
    size(B)[1] == size(B)[2] == n || throw(ArgumentError("B has wrong size"))    
    size(C)[1] == size(C)[2] == m || throw(ArgumentError("C has wrong size"))
    size(D, 1) == size(M, 1) && size(D, 2) == size(M, 2) || throw(ArgumentError("D must be the same size as M"))    

    LS_M, LS_b = system_no_null(M, M_id, B, C, D)
    LS_M, = pinv_qr(LS_M) 

    LS_x = LS_M * LS_b

    Y = unknown_to_matrix(LS_x, M_id)
    return Y
end


function system_no_null(M::AbstractArray{T}, M_id::AbstractSparseMatrixCSC,
    B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix) where{T}
    m = size(M, 1)
    n = size(M, 2)

    size(M_id, 1) == size(M, 1) && size(M_id, 2) == size(M, 2) || 
    throw(ArgumentError("M_id must be the same size as M"))
    0 <= minimum(M_id) || throw(DomainError(minimum(M_id), "M_id contains negative integers"))
    size(B)[1] == size(B)[2] == n || throw(ArgumentError("B has wrong size"))    
    size(C)[1] == size(C)[2] == m || throw(ArgumentError("C has wrong size"))
    size(D, 1) == size(M, 1) && size(D, 2) == size(M, 2) || throw(ArgumentError("D must be the same size as M"))     

    N = length(unique(M_id.nzval))

    LS_A = zeros(T, (N, N))
    LS_b = zeros(T, N)

    for idx in findall(!iszero, M_id)
        J, L = findnz(M_id[Tuple(idx)[1], :])
        k = 1
        for l in L
            @views LS_A[M_id[idx], l] = LS_A[M_id[idx], l] + B[J[k], Tuple(idx)[2]]
            k += 1
        end

        J, L = findnz(M_id[:, Tuple(idx)[2]])
        k = 1
        for l in L
            @views LS_A[M_id[idx], l] = LS_A[M_id[idx], l] + C[Tuple(idx)[1], J[k]]
            k += 1
        end
    end

    @views for i in 1:N
        LS_b[i] = sum(D[M_id .== i])
    end

    return LS_A, LS_b
end


function system_no_null(M::AbstractArray{<:Int}, M_id::AbstractSparseMatrixCSC, 
    pinv_MTM::AbstractMatrix, pinv_MMT::AbstractMatrix)
    return system_no_null(float(M), M_id, pinv_MTM, pinv_MMT)
end


function unknown_to_matrix(x::AbstractVector, M_id::AbstractSparseMatrixCSC)
    m = size(M_id, 1)
    n = size(M_id, 2)

    M_id_i, M_id_j, M_id_v = findnz(M_id)
    return sparse(M_id_i, M_id_j, x[M_id_v], m, n)
end


function sps_impose_action!(Y::AbstractSparseMatrixCSC, M::AbstractArray{T}, R_mat::AbstractMatrix, L_mat::AbstractMatrix,
     tol=eps(T)::Real, max_iters=300::Integer) where{T}
    tol > 0 || throw(DomainError("tol has to be greater than zero"))
    size(M, 1) == size(Y, 1) && size(M, 2) == size(Y, 2) || throw(ArgumentError("Y must be the same size as"))

    m = size(Y, 1)
    n = size(Y, 2)
    nR1, nR2 = size(R_mat)
    nL1, nL2 = size(L_mat)

    n == nR1 || throw(ArgumentError("Y * R_mat not possible"))
    m == nL1 || throw(ArgumentError("Y' * L_mat not possible"))

    # Quantities don't change when iterating
    norm_R_mat = norm(R_mat)
    norm_L_mat = norm(L_mat)

    M_R_mat = M * R_mat
    M_L_mat = M' * L_mat

    nnz_Y = count(!iszero, Y)
    Y_i, Y_j, _ = findnz(Y)

    # Quantities that do change when iterating
    Lag_R = zeros(T, m, nR2)
    Lag_L = zeros(T, n, nL2)

    d_R = Y * R_mat - M_R_mat
    d_L = Y' * L_mat - M_L_mat

    q_R = - d_R
    q_L = - d_L

    d_R_proj_vec = zeros(T, nnz_Y, 1)
    d_L_proj_vec = zeros(T, nnz_Y, 1)
    d_proj_vec = zeros(T, nnz_Y, 1)

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


function sps_impose_action!(Y::AbstractSparseMatrixCSC, M::AbstractArray{<:Int}, R_mat::AbstractMatrix, 
    L_mat::AbstractMatrix, tol::Real, max_iters=1000::Integer)
    return sps_impose_action!(Y, float(M), R_mat, L_mat, eps(eltype(M)), max_iters)
end

