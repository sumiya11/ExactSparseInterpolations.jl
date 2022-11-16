#=
    Newton algorithm.

    Univariate, Dense, Deterministic, Blackbox.

    Can be turned into an adaptive algorithm.

    d+1 evaluation points, where d is the polyomial degree.
=#

mutable struct Newton{Ring, Poly} <: AbstractPolynomialInterpolator
    ring::Ring
    d::Int
    f::Poly
    q::Poly
    i::Int
    function Newton(ring::Ring; d::Union{Integer, Nothing}=nothing) where {Ring}
        @assert Nemo.nvars(ring) == 1
        @assert isnothing(d) || d >= 0
        isnothing(d) && (d = typemax(Int))
        new{Ring, elem_type(ring)}(ring, d, zero(ring), zero(ring), 0)
    end
end

function Base.copy(n::Newton)
    Newton(n.ring, n.d, n.f, n.q, n.i)
end

function next_point!(N::Newton)
    error("next! should be used directly.")
end

function next!(n::Newton, x, y)
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
        m = y - n.f(x)
        n.i > n.d && @assert iszero(m)
        n.f = n.f + inv(n.q(x))*n.q*m
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
