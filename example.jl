# using ExactSparseInterpolations
using ExactInterpolations
using Nemo

# declare Q[x,y,z]
R, (x, y, z) = QQ["x","y","z"]

# example rational function
interpolate_me = Blackbox((x^2 - 4y*z + 3)//(x*y^2 + z))

# interpolator over R with numerator and denominator degrees
# not exceeding 3 and 4 respectively
vdhl = vanDerHoevenLecerf(R, 3, 4)

i = interpolate!(vdhl, f)

@assert i == interpolate_me
