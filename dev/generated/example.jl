using SubspacePreservingSparsification
using SparseArrays
using LinearAlgebra
M = [17.05 16.98 0.3 6.99 7; 16.98 0.2 7.1 6.9 0; 0.3 7.1 -12 0.01 17; 6.99 6.9 0.01 -11.97 0; 7 0 17 0 -0.1]

M_id = sparsity_pattern(M, 0.6, 2)

bin_sparse_matrix!(M, M_id, 200)

binned_minimization(M_id, zeros(5, 5), Matrix{Float64}(I, 5, 5), M)

sparsify(M, 0.6, 2, 200)

sparsify(M, 0.6, 2, 200, true)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

