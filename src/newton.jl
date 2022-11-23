#=
    Newton polynomial interpolation algorithm.

    Univariate, Dense, Deterministic, Blackbox.

    Can be turned into an adaptive algorithm.

    d+1 evaluation points, where d is the polyomial degree.
    O(d^2logd) arithmetic operations in the ground field.
=#

mutable struct Newton{Ring, Poly} <: AbstractPolynomialInterpolator
    ring::Ring
    d::Int
    f::Poly
    q::Poly
    i::Int
    function Newton{Ring,Poly}(ring::Ring, d::Union{Integer, Nothing}) where {Ring, Poly}
        @assert Nemo.nvars(ring) == 1
        @assert isnothing(d) || d >= 0
        isnothing(d) && (d = typemax(Int))
        new{Ring, elem_type(ring)}(ring, d, zero(ring), zero(ring), 0)
    end
end

Newton(ring::Ring; d::Union{Integer, Nothing}=nothing) where {Ring} = Newton{Ring, elem_type(ring)}(ring, d)

function Base.copy(n::Newton{Ring,Poly}) where {Ring,Poly}
    Newton{Ring,Poly}(n.ring, n.d, n.f, n.q, n.i)
end

function Base.empty!(n::Newton)
    n.i = 0
    n.q = zero(n.ring)
    n.f = zero(n.ring)
end

# Performs a step in the newton interpolator,
# given the new interpolation point: x and y = f(x)
function next!(n::Newton, x::T, y::T) where {T}
    R = n.ring
    z = gen(R)
    success = false
    if iszero(n.i)
        # if handling first point in the sequence
        n.f = R(y)
        n.q = z - x
        n.i = 1
        success = false
    else
        # O(dlogd)
        m = y - n.f(x)
        # O(dlogd)
        n.f = n.f + inv(n.q(x))*n.q*m
        # O(d)
        n.q = n.q*(z - x)
        n.i += 1
        success = iszero(m) || n.i > n.d
    end
    success, n.f
end

function interpolate!(n::Newton, blackbox)
    K = base_ring(n.ring)
    f = one(n.ring)
    success = false
    i = 0
    while i <= n.d && !success
        x = random_point(K)
        y = blackbox(x)
        success, f = next!(n, x, y)
        i += 1
    end
    f
end

function interpolate!(n::Newton, xs::Vector{T}, ys::Vector{T}) where {T}
    f = one(n.ring)
    i = 1
    while i <= length(xs)
        x = xs[i]
        y = ys[i]
        _, f = next!(n, x, y)
        i += 1
    end
    f
end
