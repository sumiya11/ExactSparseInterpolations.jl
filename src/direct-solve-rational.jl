#=
    Direct dense univariate rational interpolation.

    Solves the system F*Q - P = 0 to find numerator P and denominator Q,
    assuming degrees of P and Q do not exceed N and D respectively.

    N + D + 1 blackbox evaluations.

    (!) Assumes the constant term in Q is nonzero, and normalizes by it.
=#

mutable struct DirectSolveRational{Ring} <: AbstractRationalInterpolator
    ring::Ring
    N::Int
    D::Int
end

function DirectSolveRational(ring::Ring, N::Integer, D::Integer) where {Ring}
    @assert N >= 0 && D >= 0 
    DirectSolveRational(ring, N, D)
end

function interpolate!(dsr::DirectSolveRational, blackbox)
    K = base_ring(dsr.ring)
    xs = map(_ -> random_point(K), 0:dsr.N + dsr.D)
    ys = map(blackbox, xs)
    interpolate!(dsr, xs, ys)
end

function interpolate!(dsr::DirectSolveRational, xs::Vector{T}, ys::Vector{T}) where {T}
    @assert length(xs) == length(ys) == dsr.N + dsr.D + 1
    R = dsr.ring
    K = base_ring(R)
    N, D = dsr.N, dsr.D
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
