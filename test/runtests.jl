using Test
using TestSetExtensions

using Nemo
import Primes
import Random

include("../src/ExactSparseInterpolations.jl")

@testset "All tests" verbose = true begin
    @includetests ["generic", "discrete-log"]
    @includetests ["fastgcd"]
    @includetests ["div-and-conq"]

    @includetests ["cauchy"]
    @includetests ["ben-or-tiwari.jl"]
end
