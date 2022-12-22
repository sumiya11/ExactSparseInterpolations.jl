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
@time num, den = interpolate!(vdhl, bb)
## prints
## 0.067879 seconds (128.57 k allocations: 8.150 MiB)

@assert num//den == func
