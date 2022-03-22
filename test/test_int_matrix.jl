using SparseArrays

A = [2 3; 1 0; 0 6]

@test @inferred sparsify(A, 0.4, Inf, 200, true) ≈ sparse([2.0 3.0; 1.0 0; 0 6.0])
@test @inferred sparsify(A', 0.4, Inf, 200, true) ≈ sparse([2.0 3.0; 1.0 0; 0 6.0]')
A_id = @inferred sparsity_pattern(A, 0.4, Inf)
A_id_t = @inferred sparsity_pattern(A', 0.4, Inf)
@test A_id ≈ sparse([1 1; 1 0; 0 1])
@test A_id_t ≈ sparse([1 1; 1 0; 0 1]')
@test @inferred bin_sparse_matrix!(A, A_id, 200) ≈ sparse([1 3; 2 0; 0 4])
@test @inferred bin_sparse_matrix!(A', A_id_t, 200) ≈ sparse([1 3 0; 2 0 4])
