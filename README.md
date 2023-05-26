# ExactSparseInterpolations.jl

A proof-of-concept implementation of new sparse polynomial factorization algorithms in Julia. 

### Installation

In Julia, run the following:

```julia
import Pkg; Pkg.dev(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
Pkg.add("Nemo")
```

### Example

This package can be used in combination with Nemo.jl.
To factor some polynomials, first, load the packages:

```julia
using ExactSparseInterpolations
using Nemo
```

Create a polynomial ring in three variables and a polynomial:

```julia
# declare Z[x,y,z]/<3*2^30+1>
R, (x, y, z) = GF(3*2^30+1)["x","y","z"]

f = (x + y^2 + z^2 + 7) * (x^2 + x*y*z + 1) * (y^2 + z)
```

Then, you can use the `top_level_factorize` function:

```julia
top_level_factorize(f)
```

As a result, you will obtain:

```julia
3-element Vector{gfp_mpoly}:
 x^2 + x*y*z + 1
 x + y^2 + z^2 + 7
 y^2 + z
```

