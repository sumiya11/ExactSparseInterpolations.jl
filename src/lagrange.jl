#=
    Interpolation via Lagrange basis polynomials.

    Univariate, Dense, Deterministic, No-Early-Termination,
    Blackbox & Given points.

    O()
=#

mutable struct Lagrange{T1} <: AbstractInterpolator
    ring::T1
    d::Int
    function Lagrange(ring::T; d::Integer=1) where {T}
        @assert Nemo.nvars(ring) == 1
        @assert d >= 0
        new{T}(ring, d)
    end
end

function Nemo.interpolate(L::Lagrange{T1}, blackbox) where {T1}
    xs = map(_ -> first(generic_point(parent(L))), 0:L.d)
    ys = map(blackbox, xs)
    interpolate(L, xs, ys)
end

function Nemo.interpolate(L::Lagrange{T1}, 
                    p::Vector{T}, m::Vector{T}) where {T1, T}
    R, x = parent(L), first(gens(L))
    (length(m) == 1) && return R(first(m))
    fx = zero(R)
    num = prod((x - pj) for pj in p) 
    for i in eachindex(p, m)
        den = prod((p[i] - p[j]) for j in eachindex(p) if j != i)
        fx += m[i]*divexact(divexact(num, x - p[i]), den)
    end
    fx
end
