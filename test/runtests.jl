using Test
using TestSetExtensions

using Nemo
import Primes
import Random

include("../src/ExactSparseInterpolations.jl")

@testset "All tests" verbose=true begin
    @includetests ["generic", "discrete-log"]
    @includetests ["fastgcd"]
    @includetests ["div-and-conq"]
    @includetests ["faster-cauchy"]
    @includetests ["faster-ben-or-tiwari"]
    @includetests ["faster-van-der-hoeven-lecerf"]

    @includetests ["hensel"]
    @includetests ["bivariate-factor-ff"]
end
