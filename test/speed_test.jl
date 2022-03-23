using SubspacePreservingSparsification
using SubspacePreservingSparsification: near_zero_row_col
using BenchmarkTools
using MAT
using SparseArrays
using DelimitedFiles
using LinearAlgebra


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

# print("Computation\n")
# @benchmark sparsify($A, 0.4, 2, max_num_bin, true)
