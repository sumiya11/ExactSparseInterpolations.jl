using Test
using TestSetExtensions

using Nemo

include("../src/ExactSparseInterpolations.jl")

@testset "All tests" verbose=true begin
    @includetests ["berlekamp-massey"]
    @includetests ["newton"]
    @includetests ["ben-or-tiwari"]
    @includetests ["cuyt-lee"]
    @includetests ["van-der-hoeven-lecerf"]
end
