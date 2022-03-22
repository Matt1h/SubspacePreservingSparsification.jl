"""
    pinv_M, rnull_optional, lnull_optional = pinv_qr(M::AbstractArray)

Computes Moore-Penrose pseudoinverse of `M` using a rank-revealing QR
factorization and optionally returns an orthonormal basis for right
and left null-space. This can be significantly faster than using
SVD for computing the pseudoinverse.  Moreover, `A` can be sparse or
dense.

# Examples
```julia-repl
julia> pinv_qr([3 8 17; 2 1 4; 8 3 21])
([-0.10344827586206914 1.344827586206895 -0.17241379310344793;0.1149425287356321 0.8390804597701144 -0.2528735632183907;0.022988505747126506 -0.6321839080459762 0.1494252873563217], Matrix{Float64}(undef, 3, 0), Matrix{Float64}(undef, 3, 0))
```
"""
function pinv_qr(M::AbstractArray{T}) where {T}
    m = size(M, 1)
    n = size(M, 2)
    
    F = qr(M)
    Q = F.Q*I  # so Q is allways m x m, with Matrix(F.Q) we get a m x n Matrix for m >= n
    R = Matrix(F.R)

    # compute rank 
    abs_diag_R = abs.(diag(R))
    tol = 100 * minimum(size(M)) * maximum(abs_diag_R) * eps(T)  # MAGIC CONSTANT
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
        rnull = zeros(T, n, 0)
    end

    # left null space
    if left_null_size > 0
        lnull = Q[: ,end - (left_null_size - 1): end]
    else
        lnull = zeros(T, m, 0)
    end

    return M_pinv, rnull, lnull
end


function pinv_qr(M::AbstractArray{<:Int})
    pinv_qr(float(M))
end