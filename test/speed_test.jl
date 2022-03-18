using SubspacePreservingSparsification
using SubspacePreservingSparsification: near_zero_row_col
using BenchmarkTools
using MAT
using SparseArrays
using DelimitedFiles
using LinearAlgebra
using ProfileView


function init_script()
    ROOT_DIR = realpath(dirname(@__FILE__) * "/..")
    cd(ROOT_DIR)
end

init_script()
sq = "square100"
file = matopen(joinpath("test", "test_data", sq, "untransformed_matrices", "centrosymmetric.mat"))
A = read(file, "A")
close(file)


max_num_bin = 50

pinv_A, rnull, lnull = pinv_qr(A)

# sparsity pattern
num_near_zero_rows, num_near_zero_cols = near_zero_row_col(A)
min_per_row = max(0, min(size(rnull)[2] - num_near_zero_cols, size(A)[2]))
min_per_col = max(0, min(size(lnull)[2] - num_near_zero_rows, size(A)[1]))

A_pat = p_norm_sparsity_matrix(A, 0.4, 2, min_per_row, min_per_col)
A_id = bin_sparse_matrix!(A, A_pat, max_num_bin)

sps_compute([1;2;0], 0.6, 2, 200)

# pinv_ATA = pinv_A * pinv_A'
# pinv_AAT = pinv_A' * pinv_A

# n = Int(size(A)[1]/2)
# z = spzeros(n, n)
# i = sparse(I, n, n)
# K = [z i; -i z]
# norm(K*X - X'*K)

# n = size(A)[1]
# z = spzeros(n-1)
# i = sparse(I, n-1, n-1)
# C = [z i; -1 z']
# norm(C*X - X*C)

# J = reverse(sparse(I, size(X)[1], size(X)[1]), dims=1)
# norm(X*J - J*X)

# print("Computation\n")
# @benchmark ssa_compute($A, 0.4, 2, max_num_bin, true)

# function speedtest()
#     for i in 1:100
#         ssa_compute(A, 0.4, 2, max_num_bin, true)
#     end
# end
# ProfileView.@profview speedtest()
# ProfileView.@profview speedtest()
