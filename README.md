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

# declare Z[x,y,z]/<3*2^30+1>
R, (x, y, z) = GF(3*2^30+1)["x","y","z"]

# create an example function for interpolation
func = (x^10 - y + z + 8)^5 // (x + 2x*y)

# obtain information about degrees/terms/partial degrees of func
info = ExactSparseInterpolations.getboundsinfo(func)

# wrap func into a blackbox
bb = ExactSparseInterpolations.Blackbox(func)

# create an interpolator, using the information of func
vdhl = ExactSparseInterpolations.FasterVanDerHoevenLecerf(R, info)

# interpolate (note that this mutates vdhl object)
@time num, den = interpolate!(vdhl, f)
## prints
## 0.101146 seconds (131.54 k allocations: 8.196 MiB)

@assert num//den == func

```
