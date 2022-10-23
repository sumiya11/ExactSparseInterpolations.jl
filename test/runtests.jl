using Test
using TestSetExtensions

using Nemo

include("../src/ExactInterpolations.jl")

@testset "All tests" begin
    @includetests ["berlekamp-massey"]
    @includetests ["lagrange","newton"]
    @includetests ["zippel"]
end
