module SubspacePreservingSparsification

using LinearAlgebra
using SparseArrays
using SparseArrays: AbstractSparseMatrixCSC

include("math_functions.jl")
include("sparse_pattern.jl")
include("bin_pattern.jl")
include("compute.jl")

export pinv_qr
export sparsity_pattern
export bin_sparse_matrix!
export binned_minimization
export sparsify

end  # end of module
