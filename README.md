# ExactSparseInterpolations.jl

A proof-of-concept implementation of some new sparse polynomial factorization algorithms in Julia. 

### Installation

In Julia, run the following:

```julia
import Pkg; Pkg.develop(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
Pkg.add("Nemo")
```

**Disclaimer.**
Recently, there was a number of breaking changes in Nemo.jl. 
The examples below should work out of the box, but if they don't, you might want to install a specific verion of Nemo.jl.
To do that, assuming you are in the root directory of this project, simply do

```julia
import Pkg; 
Pkg.activate("version_freeze"); 
Pkg.instantiate()
```

### Example

This package can be used in combination with Nemo.jl.
First, load the packages:

```julia
using ExactSparseInterpolations
using Nemo
```

Then, you can use the function `irreducible_factorization` over a finite field:

```julia
# declare Z[x,y,z]/<3*2^30+1>
R, (x1, x2, x3) = GF(3*2^30+1)["x1","x2","x3"]

F = (x1 + x2)*(x2 + x3)*(x1*x2 - x1*x3 - 2)*(x1*x3 - x2*x3 - 3)*(x1*x2*x3 - 99)

c, factors = irreducible_factorization(F)
```

As a result, you will obtain:

```julia
(1, gfp_mpoly[x1*x3 + 3221225472*x2*x3 + 3221225470, x1*x2 + 3221225472*x1*x3 + 3221225471, x1*x2*x3 + 3221225374, x2 + x3, x1 + x2])
```

You can also use the version with a tri-variate substitution:

```julia
c_2, factors_2 = irreducible_factorization_2(F)

@assert c_2 == c && sort(factors_2, by=leading_monomial) == sort(factors, by=leading_monomial)
```