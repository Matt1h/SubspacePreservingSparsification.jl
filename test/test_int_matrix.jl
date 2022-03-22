A = [2 3; 1 0; 0 6]

@test sparsify(A, 0.4, Inf, 200, true) ≈ sparse([2.0 3.0; 1.0 0; 0 6.0])