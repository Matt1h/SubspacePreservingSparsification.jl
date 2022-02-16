using Revise
using SSA
using BenchmarkTools
using MAT
using SparseArrays
using DelimitedFiles
using LinearAlgebra
using ProfileView

file = matopen(joinpath("test", "test_data", "square100", 
"untransformed_matrices", "with_nullspace.mat"))
A = read(file, "A")
close(file)

max_num_bin = 50

pinv_A, rnull, lnull = pinv_qr(A)

# sparsity pattern
num_near_zero_rows, num_near_zero_cols = SSA.near_zero_row_col(A)
min_per_row = max(0, min(size(rnull)[2] - num_near_zero_cols, size(A)[2]))
min_per_col = max(0, min(size(lnull)[2] - num_near_zero_rows, size(A)[1]))

A_pat = p_norm_sparsity_matrix(A, 0.4, 2, min_per_row, min_per_col)
A_id = bin_sparse_matrix!(A, A_pat, max_num_bin)

pinv_ATA = pinv_A * pinv_A'
pinv_AAT = pinv_A' * pinv_A

LS_A, LS_b = SSA.ssa_system_no_null(A, A_id, pinv_ATA, pinv_AAT)

# LS_A, = pinv_qr(LS_A)
# LS_x = LS_A * LS_b

# print("Sparsity Pattern\n")
# @benchmark p_norm_sparsity_matrix(A, 0.4, 2, min_per_row, min_per_col)

# print("Binning Pattern\n")
# @benchmark bin_sparse_matrix(A, A_pat, 0)

# print("Computation\n")
# @benchmark ssa_compute(A, 0.4, 2, max_num_bin, true)

# ssa_compute(A, 0.4, 2, max_num_bin, true)
# @time ssa_compute(A, 0.4, 2, max_num_bin, true)

# print("pseudo inverse\n")
# @benchmark pinv_qr(LS_A)

# print("minimization\n")
# @benchmark SSA.ssa_minimization(A, A_id, pinv_A)

# print("find equation system\n")
# @benchmark SSA.ssa_system_no_null(A, A_id, pinv_ATA, pinv_AAT)

# print("unknown to matrix\n")
# @benchmark SSA.ssa_unknown_to_matrix(LS_x, A_id)

# @benchmark pinv_qr(LS_A)

# X = SSA.ssa_minimization(A, A_id, pinv_A)
# @benchmark SSA.ssa_impose_action!(X, A, rnull, lnull)


# function speedtest()
#     for i in 1:100
#     ssa_compute(A, 0.4, 2, max_num_bin, true)
#     end
# end
# ProfileView.@profview speedtest()
# ProfileView.@profview speedtest()


# @code_warntype ssa_compute(A, 0.4, 2, max_num_bin, true)

# A =     [16     2     3    13 1
# 5    11    10     8 2 
# 9     7     6    12 0 
# 4    14    15     1 0]'
# print(A)
# print(nullspace(A))

# function default_ssa_compute_bin_no_null(M)
#     M = ssa_compute(M, 0.4, 2, 50)
#     return M
# end


# function default_ssa_compute_no_bin_no_null(M)
#     M = ssa_compute(M, 0.4, 2, 0)
#     return M
# end


# function my_default_ssa_compute_bin_null(M)
#     M = ssa_compute(M, 0.4, 2, 50, true)
#     return M
# end


# function my_default_compute_sparse_pattern(M)
#     pinv_M, rnull, lnull = pinv_qr(M)

#     # sparsity pattern
#     num_near_zero_rows, num_near_zero_cols = SSA.near_zero_row_col(M)
#     min_per_row = max(0, min(size(rnull)[2] - num_near_zero_cols, size(M)[2]))
#     min_per_col = max(0, min(size(lnull)[2] - num_near_zero_rows, size(M)[1]))
#     M = p_norm_sparsity_matrix(M, 0.4, 2, min_per_row, min_per_col)
#     return M
# end


# function my_default_compute_bin_pattern(M)
#     pinv_M, rnull, lnull = pinv_qr(M)

#     # sparsity pattern
#     num_near_zero_rows, num_near_zero_cols = SSA.near_zero_row_col(M)
#     min_per_row = max(0, min(size(rnull)[2] - num_near_zero_cols, size(M)[2]))
#     min_per_col = max(0, min(size(lnull)[2] - num_near_zero_rows, size(M)[1]))
#     M_pat = p_norm_sparsity_matrix(M, 0.4, 2, min_per_row, min_per_col)
#     M = bin_sparse_matrix(M, M_pat, 50)
#     return M
# end


# function init_script()
#     ROOT_DIR = realpath(dirname(@__FILE__) * "/..")
#     cd(ROOT_DIR)
# end

# init_script()

# type_names = readdlm(joinpath("test", "test_data", "square100", "type_names.txt"))
# trans_names = readdlm(joinpath("test", "test_data", "square100", "trans_names.txt"))


# transformations = Dict(
#     "pinv_rrqr" => pinv_qr,
#     "default_ssa_compute_bin_no_null" => default_ssa_compute_bin_no_null,
#     "default_ssa_compute_no_bin_no_null" => default_ssa_compute_no_bin_no_null,
#     "default_ssa_compute_bin_null" => my_default_ssa_compute_bin_null,
#     "sparse_pattern" => my_default_compute_sparse_pattern,
#     "bin_pattern" => my_default_compute_bin_pattern,
# )


# file = matopen(joinpath("test", "test_data", "square10", "untransformed_matrices", "with_nullspace.mat"))
# A = read(file, "A")
# close(file)
# for j_trans_name in trans_names
#     trans = transformations[j_trans_name]
#     trans_A = trans(A)
#     file = matopen(joinpath("test", "test_data", "square10", 
#     "transformed_matrices", "matlab", "with_nullspace", j_trans_name*".mat"))
#     trans_A_mat = read(file, "trans_A")  # TODO: no sparse matrix
#     close(file)
# end
