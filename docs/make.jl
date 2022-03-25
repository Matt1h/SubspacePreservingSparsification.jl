using Documenter
using SubspacePreservingSparsification
using Literate

EXAMPLE_DIR = joinpath(@__DIR__, "literate")
OUT_DIR = joinpath(@__DIR__, "src\\generated")

# Use Literate.jl to generate docs and notebooks of examples
for example in readdir(EXAMPLE_DIR)
    EXAMPLE = joinpath(EXAMPLE_DIR, example)

    Literate.markdown(EXAMPLE, OUT_DIR; documenter=true) # markdown for Documenter.jl
    Literate.notebook(EXAMPLE, OUT_DIR) # .ipynb notebook
    Literate.script(EXAMPLE, OUT_DIR) # .jl script
end


makedocs(;
    modules=[SubspacePreservingSparsification],
    sitename="SubspacePreservingSparsification.jl",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true", assets=String[]),
    pages=[
        "Home" => "index.md",
        "example" => "generated/example.md",
    ],
)

deploydocs(; repo="github.com/Matt1h/SubspacePreservingSparsification.jl")