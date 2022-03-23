using SubspacePreservingSparsification
using Test

@testset "SubspacePreservingSparsification.jl" begin
    @testset "Bin pattern" begin
        include("test_bin_pattern.jl")
    end
    @testset "diffenrent types without nullspace" begin
        include("test_compare_with_matlab.jl")
    end
    @testset "with nullspace" begin
        include("test_compare_with_matlab_nullspace.jl")
    end
    @testset "Integer Matrix" begin
        include("test_int_matrix.jl")
    end
end
