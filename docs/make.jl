push!(LOAD_PATH,"../src/")

using Documenter
using SSA

makedocs(modules=[SSA],
sitename="My Documentation")