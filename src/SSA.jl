module SSA

using StatsBase
using SparseArrays
using LinearAlgebra

include("math_functions.jl")
include("sparse_pattern.jl")
include("bin_pattern.jl")
include("compute.jl")

export pinv_qr
export p_norm_sparsity_matrix
export bin_sparse_matrix!
export ssa_compute

end  # end of module
