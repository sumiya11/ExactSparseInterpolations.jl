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