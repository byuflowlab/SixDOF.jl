using SixDOF
using Test

@testset "Forward/Reverse Diff" begin
    include("testderiv.jl")
end

@testset "Reference frames" begin
    include("testframes.jl")
end
