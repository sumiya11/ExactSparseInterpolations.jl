# ExactSparseInterpolations.jl

A proof-of-concept implementation of some new sparse polynomial factorization algorithms in Julia. 

### Installation

In Julia, run the following:

```julia
import Pkg; Pkg.develop(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
Pkg.add("Nemo")
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
[ Info: Using random dilation gfp_mpoly[x1, x2, x3] -> gfp_mpoly[613987896*x1, 2209294782*x2, 1503459614*x3]
┌ Info: Input polynomial is in 3 variables over Galois field with characteristic 3221225473.
│ Performing the substitution xᵢ -> ξʷⁱ xᵢ ((1 - λ) + αᵢλ) (t + βᵢ u + γᵢ v))
│ Using the following evaluation points:
│ α = gfp_elem[2, 3, 5]
│ β = gfp_elem[3073480812, 26366625, 1138354345]
│ γ = gfp_elem[48733866, 2612364570, 1542073172]
│ And regularizing weight:
└ w = [4, 1, 4]
[ Info: Computing the factorizations at indices 2..2
[ Info: Trying to interpolate from 2 evaluations at powers of α = gfp_elem[2, 3, 5]
[ Info: Computing the factorizations at indices 3..4
[ Info: Trying to interpolate from 4 evaluations at powers of α = gfp_elem[2, 3, 5]
[ Info: Successful interpolation!
(1, gfp_mpoly[x1*x3 + 3221225472*x2*x3 + 3221225470, x1*x2 + 3221225472*x1*x3 + 3221225471, x1*x2*x3 + 3221225374, x2 + x3, x1 + x2])
```
