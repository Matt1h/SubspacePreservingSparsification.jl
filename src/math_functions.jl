function pinv_qr(M::AbstractArray{T}) where {T}
    m = size(M, 1)
    n = size(M, 2)
    
    F = qr(M)
    Q = F.Q*Matrix(I ,m, m)  # so Q is allways m x m, with Matrix(F.Q) we get a m x n Matrix for m >= n
    R = Matrix(F.R)

    # compute rank 
    abs_diag_R = abs.(diag(R))
    tol = 100 * minimum(size(M)) * maximum(abs_diag_R) * eps(T);  # MAGIC CONSTANT
    rnk = sum(tol .< abs_diag_R)

    # compute pseudoinverse
    M_pinv = Matrix{T}(undef, n, m)
    S = R[1:rnk, 1:rnk] \ R[1:rnk, rnk+1:end]
    W = R[1:rnk, 1:rnk] \ Q[:, 1:rnk]'
    X = (S' * S + I) \ (S' * W) 
    M_pinv[rnk+1:end, :] = X
    M_pinv[1:rnk, :] = W - S*X

    # computing null spaces
    right_null_size = n - rnk
    left_null_size = m - rnk

    # right null space 
    if right_null_size > 0
        rnull, = qr(Matrix([S; -I(right_null_size)]))
        rnull = Matrix(rnull)
    else
        rnull = zeros(n, 0)
    end

    # left null space
    if left_null_size > 0
        lnull = Q[: ,end - (left_null_size - 1): end]
    else
        lnull = zeros(m, 0)
    end

    return M_pinv, rnull, lnull
end