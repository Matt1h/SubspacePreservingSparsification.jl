```@meta
CurrentModule = SubspacePreservingSparsification
end
```
# SubspacePreservingSparsification.jl
Implementation of an algorithm, that takes a real Matrix M and finds a sparse approximation of the same size. The algorithm was developed by Chetan Jhurani under the name sparse spectral approximation (SSA). See [https://github.com/cjhurani/txssa](https://github.com/cjhurani/txssa) for more detailed documentation, also with regard to the mathematical background, and for implementations in C/C++ and in Matlab.

## Installation
To install this package and its dependencies, open the Julia REPL and run 
```julia-repl
julia> ]add SubspacePreservingSparsification
```

## Manual
```@contents
Pages = [
    "generated/example.md",
]
Depth = 1
```

## API reference
```@autodocs
Modules = [SubspacePreservingSparsification]
Order   = [:function, :type]
```