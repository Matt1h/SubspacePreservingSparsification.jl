using Documenter
using SubspacePreservingSparsification

DocMeta.setdocmeta!(
    SubspacePreservingSparsification,
    :DocTestSetup,
    :(using SparseArrays, SubspacePreservingSparsification);
    recursive=true
)

makedocs(
    modules=[SubspacePreservingSparsification],
    sitename="SubspacePreservingSparsification.jl",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true", assets=String[]),
    pages=["Getting started" => "index.md"],
    strict=:doctest,
)
deploydocs(; repo="github.com/Matt1h/SubspacePreservingSparsification.jl")
