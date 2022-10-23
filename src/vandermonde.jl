#=
    Interpolation via Vandermonde matrix inversion.

    Two algorithms:

    Univariate, Dense, Deterministic, No-Early-Termination,
    Blackbox & Given points.

    O()
=#

mutable struct Vandermonde{T1} <: AbstractInterpolator
    ring::T1
    d::Int
    function Vandermonde(ring::T; d::Integer=1) where {T}
        @assert Nemo.nvars(ring) == 1
        @assert d >= 0
        new{T}(ring, d)
    end
end

function Nemo.interpolate(V::Vandermonde{T1}, blackbox) where {T1}
    xs = map(_ -> first(generic_point(parent(L))), 0:L.d)
    ys = map(blackbox, xs)
    interpolate(L, xs, ys)
end

function vandermonde_inverse(V)

end