# ExactSparseInterpolations.jl

### Installation

Run in Julia:

```julia
import Pkg; Pkg.add(url="https://github.com/sumiya11/ExactSparseInterpolations.jl")
```

### How to use ExactSparseInterpolations.jl ?

See file `example.jl` or the simple example below:

```julia
using ExactSparseInterpolations
using Nemo

# declare Q[x,y,z]
R, (x, y, z) = QQ["x","y","z"]

# example rational function
f = Blackbox((x^2 - 4y*z + 3)//(x*y^2 + z))

# van Der Hoeven and Lecerf algorithm in R with numerator
# and denominator degrees not exceeding 3 and 4 respectively
vdhl = vanDerHoevenLecerf(R, 3, 4)

i = interpolate!(vdhl, f)
(x^2 - 4*y*z + 3, x*y^2 + z)
```

Note that the `vdhl` object is mutated in `interpolate!`.

Or, over a finite field:

```julia
# declare Z[x,y]/<2^31-1>
R, (x, y) = GF(2^31-1)["x","y"]

f = Blackbox((x - y + 8)^5 // (x + 2x*y))

vdhl = vanDerHoevenLecerf(R, 5, 2)

num, den = interpolate!(vdhl, f)

@assert num//den == (x - y + 8)^5 // (x + 2x*y)

```

Note that this will uncontrollably produce wrong results in the case the field characteristic $p$ is too small, $p$ should be at least $\Omega(n^d)$.
