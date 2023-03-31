# All interpolators are a subtype of AbstractInterpolator
abstract type AbstractInterpolator end

Nemo.parent(x::AbstractInterpolator) = getfield(x, :ring)
Nemo.base_ring(x::AbstractInterpolator) = base_ring(parent(x))
Nemo.gen(x::AbstractInterpolator) = gen(parent(x))
Nemo.gens(x::AbstractInterpolator) = gens(parent(x))
Nemo.ngens(x::AbstractInterpolator) = ngens(parent(x))

abstract type AbstractPolynomialInterpolator <: AbstractInterpolator end
abstract type AbstractRationalInterpolator <: AbstractInterpolator end
