
mutable struct RationalDifferenceNormalized{T1} <: AbstractInterpolator
    ring::T1
    N::Int
    D::Int
end

function RationalDifferenceNormalized(ring, N, D)
    RationalDifferenceNormalized(ring, N, D)
end

function Nemo.interpolate(drd::RationalDifferenceNormalized{T1}, blackbox) where {T1}
    xs = map(_ -> first(generic_point(parent(drd))), 0:drd.N + drd.D)
    ys = map(blackbox, xs)
    interpolate(drd, xs, ys)
end

function Nemo.interpolate(drd::RationalDifferenceNormalized{T1}, xs, ys) where {T1}
    @assert length(xs) == length(ys) == drd.N + drd.D + 1
    R = drd.ring
    K = base_ring(R)
    N, D = drd.N, drd.D
    S = MatrixSpace(K, N + D + 1, N + D + 1)
    B = MatrixSpace(K, N + D + 1, 1)
    P = zero(S)
    b = zero(B)
    for i in 0:N + D
        ξi = xs[i+1]
        Fi = ys[i+1]
        Qs = [Fi*ξi^j for j in 1:D]
        Ps = [-ξi^j for j in 0:N]
        P[i+1, 1:D] = Qs
        P[i+1, D+1:N+D+1] = Ps
        b[i+1, 1] = -Fi
    end
    ab = collect(Nemo.solve(P, b))
    a = ab[D+1:N+D+1]
    b = append!([one(K)], ab[1:D])
    R(a), R(b)
end
