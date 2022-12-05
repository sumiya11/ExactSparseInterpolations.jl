using Test
using TestSetExtensions

using Nemo

include("../src/ExactSparseInterpolations.jl")

@testset "All tests" verbose=true begin
    @includetests ["fastgcd"]
    @includetests ["berlekamp-massey"]
    @includetests ["newton"]
    @includetests ["ben-or-tiwari"]
    @includetests ["direct-solve-rational", "cauchy"]
    @includetests ["cuyt-lee"]
    @includetests ["van-der-hoeven-lecerf"]
    @includetests ["adaptive-van-der-hoeven-lecerf"]
end
