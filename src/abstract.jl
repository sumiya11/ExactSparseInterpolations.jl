
abstract type AbstractInterpolator end

Nemo.parent(x::AbstractInterpolator) = x.ring
Nemo.base_ring(x::AbstractInterpolator) = base_ring(parent(x))
Nemo.gen(x::AbstractInterpolator) = gen(parent(x))
Nemo.gens(x::AbstractInterpolator) = gens(parent(x))
Nemo.ngens(x::AbstractInterpolator) = ngens(parent(x))

