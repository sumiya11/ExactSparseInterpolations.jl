#=
    Interpolation via Newton method.

    Univariate, Dense, Deterministic, No-Early-Termination,
    Blackbox & Given points.

    *** Can be run incrementally point by point

    O()
=#

mutable struct Newton{T1} <: AbstractInterpolator
    ring::T1
    d::Int
    function Newton(ring::T; d::Int=1) where {T}
        @assert Nemo.nvars(ring) == 1
        @assert d >= 0
        new{T}(ring, d)
    end
end

function Nemo.interpolate(N::Newton{T1}, blackbox) where {T1}
    xs = map(_ -> first(generic_point(parent(N))), 0:N.d)
    ys = map(blackbox, xs)
    interpolate(N, xs, ys)
end

function Nemo.interpolate(N::Newton{T1}, 
                    p::Vector{T}, m::Vector{T}) where {T1, T}
    k = length(m)
    R, x = parent(N), gen(N)
    fx = R(m[1])
    qx = x - p[1]
    for i in 2:k
        fx = fx + inv(qx(p[i]))*qx*(m[i] - fx(p[i]))
        qx = (x - p[i])*qx
    end
    fx
end