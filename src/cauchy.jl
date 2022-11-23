#=
    Cauchy univariate rational interpolation.
=#

# given g and f (|f| >= |g|),
# computes and returns the triple (r, t, s), such that
# r = t*g + s*f, |r| < k, |r| is the maximal possible
function constrainedEEA(g, f, k::Integer)
    @assert degree(f) >= degree(g)
    R = parent(g)
    U = (one(R), zero(R), f)
    V = (zero(R), one(R), g)
    # in Nemo, degree(0) is -1
    while degree(V[3]) >= k
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end
    (V[3], V[2], V[1])
end

mutable struct Cauchy{Ring, I1<:AbstractPolynomialInterpolator} <: AbstractRationalInterpolator
    ring::Ring
    N::Int
    D::Int
    dense_polynomial_interpolator::I1
end

function Cauchy(ring::Ring, N::Integer, D::Integer;
        dense_polynomial_interpolator=Newton(ring, d=N + D + 2)) where {Ring}
    @assert N >= 0 && D >= 0 
    Cauchy(ring, N, D, dense_polynomial_interpolator)
end

function Base.empty!(c::Cauchy)
    empty!(c.dense_polynomial_interpolator)
end

function interpolate!(c::Cauchy, blackbox)
    K = base_ring(c.ring)
    xs = map(_ -> random_point(K), 0:c.N + c.D)
    ys = map(blackbox, xs)
    interpolate!(c, xs, ys)
end

function interpolate!(c::Cauchy, xs::Vector{T}, ys::Vector{T}) where {T}
    @assert length(xs) == length(ys) == c.N + c.D + 1
    R = c.ring
    z = gen(R)
    dpi = c.dense_polynomial_interpolator
    # F is the polynomial of degree at max N + D + 1
    # that passes through intrpolation points
    F = interpolate!(dpi, xs, ys)
    modulo = prod(z - x for x in xs)
    k = c.N
    r, t, s = constrainedEEA(F, modulo, k)
    !isunit(gcd(r, t)) && throw("Cauchy interpolation fail.")
    normfactor = leading_coefficient(t)
    divexact(r, normfactor), divexact(t, normfactor)
end
