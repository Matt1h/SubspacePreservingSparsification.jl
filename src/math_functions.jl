function pinv_qr(M::AbstractMatrix)
    m,n = size(M)
    F = qr(M)

    # compute rank 
    rnk = rank(M)

    M_pinv = pinv(M)

    # computing null spaces
    right_null_size = n - rnk
    left_null_size = m - rnk

    # right null space 
    rnull = nullspace(M)

    # left null space
    lnull = F.Q[: ,end - (left_null_size - 1): end]
    return M_pinv, rnull, lnull
end