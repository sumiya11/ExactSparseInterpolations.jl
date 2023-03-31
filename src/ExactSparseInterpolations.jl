module ExactSparseInterpolations

import AbstractAlgebra
using Nemo
using Primes
import Random
import Combinatorics, Permutations
using LinearAlgebra

include("general/abstract.jl")

include("general/generic.jl")
include("general/discrete-log.jl")
include("general/field-generators.jl")
include("general/div-and-conq.jl")
include("general/blackbox.jl")
include("general/fastgcd.jl")
include("general/utils.jl")

include("interpolation/berlekamp-massey.jl")
include("interpolation/dft.jl")
include("interpolation/newton.jl")
include("interpolation/ben-or-tiwari.jl")
include("interpolation/faster-ben-or-tiwari.jl")
include("interpolation/javadi-monagan.jl")
include("interpolation/cauchy.jl")
include("interpolation/faster-cauchy.jl")
include("interpolation/direct-solve-rational.jl")
include("interpolation/cuyt-lee.jl")
include("interpolation/van-der-hoeven-lecerf.jl")
include("interpolation/faster-van-der-hoeven-lecerf.jl")
include("interpolation/adaptive-van-der-hoeven-lecerf.jl")
include("interpolation/v-d-h-l-gcd.jl")

include("factorization/univariate-factor-ff.jl")
include("factorization/hensel.jl")
include("factorization/bivariate-factor-ff.jl")
include("factorization/multivariate-factor-ff.jl")

export Blackbox

export BerlekampMassey
export Newton
export BenOrTiwari

export CuytLee

export vanDerHoevenLecerf
export vanDerHoevenLecerfGCD
export adaptiveVanDerHoevenLecerf

export interpolate!
export next!
export next_point!

end
