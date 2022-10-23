module ExactInterpolations

import AbstractAlgebra
using Nemo
using Primes

include("abstract.jl")

include("generic.jl")

include("blackbox.jl")

include("berlekamp-massey.jl")

include("vandermonde.jl")

include("lagrange.jl")

include("newton.jl")

include("zippel.jl")

include("ben-or-tiwari.jl")

include("rational-difference.jl")

include("cuyt-lee.jl")

export BerlekampMassey
export Lagrange
export Newton
export Zippel
export BenOrTiwari

export CuytLee

export interpolate

end