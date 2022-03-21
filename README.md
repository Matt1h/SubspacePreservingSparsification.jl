# SubspacePreservingSparsification.jl
Implementation of an algorithm, that takes a real Matrix M and finds a sparse approximation of the same size. The algorithm was developed by Chetan Jhurani under the name sparse spectral approximation (SSA). See [https://github.com/cjhurani/txssa](https://github.com/cjhurani/txssa) for more detailed documentation, also with regard to the mathematical background, and for implementations in C/C++ and in Matlab.

## Installation
To install this package and its dependencies, open the Julia REPL and run 
```julia-repl
julia> ]add SubspacePreservingSparsification
```

## Example

```jldoctest
julia> ssa_compute([16.99 65; 0.1 17.01], 0.6, 2, 200)
2×2 SparseMatrixCSC{Float64, Int64} with 3 stored entries:
 16.8041  64.2499
   ⋅      16.8041
```