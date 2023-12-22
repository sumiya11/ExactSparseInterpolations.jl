# ExactSparseInterpolations.jl

A proof-of-concept implementation of some new sparse polynomial factorization algorithms in Julia. 

### Installation

In Julia, run the following:

```julia
import Pkg; Pkg.dev(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
Pkg.add("Nemo")
```

### Example

This package can be used in combination with Nemo.jl.
First, load the packages:

```julia
using ExactSparseInterpolations
using Nemo
```

Then, you can use the function `irreducible_factorization`:

```julia
# declare Z[x,y,z]/<3*2^30+1>
R, (x1, x2, x3) = GF(3*2^30+1)["x1","x2","x3"]

F = (x1 + x2)*(x2 + x3)*(x1*x2 - x1*x3 - 2)*(x1*x3 - x2*x3 - 3)*(x1*x2*x3 - 99)

c, factors = irreducible_factorization(F)
```

As a result, you will obtain:

```julia

```
