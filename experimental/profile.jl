using BenchmarkTools, Random
using Nemo, Primes
# using ExactSparseInterpolations

include("../testsuite.jl")

F = GF(2^62 + 135)

nterms = 30
nv = 5
deg = 200

func = random_rational_function(1:2^16, nv, deg, deg, nterms, nterms, ground_field=F);
R = base_ring(parent(func))
metainfo = ExactSparseInterpolations.getboundsinfo(func)
blackbox = ExactSparseInterpolations.Blackbox(func);

vdhl = ExactSparseInterpolations.FasterVanDerHoevenLecerf(R, metainfo);

@profview ExactSparseInterpolations.interpolate!(vdhl, blackbox);

