#=
    Cauchy univariate rational interpolation.
=#

mutable struct Cauchy{Ring, I1<:AbstractPolynomialInterpolator} <: AbstractRationalInterpolator
    ring::Ring
    N::Int
    D::Int
    dense_polynomial_interpolator::I1
end

function Cauchy(ring::Ring, N::Integer, D::Integer; dense_polynomial_interpolator=Newton(ring, d=N + D + 2)) where {Ring}
    @assert N >= 0 && D >= 0 
    Cauchy(ring, N, D, dense_polynomial_interpolator)
end

# given (polynomials) g and f (|f| >= |g|),
# computes and returns a single row from the EEA algorithm (r, t, s), 
# such that r = t*g + s*f, |r| < k, where |r| is the maximal possible
function constrainedEEA(g, f, k::Integer)
    @assert degree(f) >= degree(g)
    R = parent(g)  # = K[x]
    U = (one(R), zero(R), f)  # = (1, 0, f)
    V = (zero(R), one(R), g)  # = (0, 1, g)
    # in Nemo, degree(0) is -1
    while degree(V[3]) > k
        q = div(U[3], V[3])
        T = U .- q .* V
        U = V
        V = T
    end
    (V[3], V[2], V[1])
end

# refreshes the current state of the interpolator
function Base.empty!(c::Cauchy)
    empty!(c.dense_polynomial_interpolator)
end

# interpolates the blackbox function using the Cauchy algorithm
function interpolate!(c::Cauchy, blackbox)
    K = base_ring(c.ring)
    xs = distinct_points(K, c.N + c.D + 2)
    ys = map(blackbox, xs)
    interpolate!(c, xs, ys)
end

function interpolate!(c::Cauchy, xs::Vector{T}, ys::Vector{T}) where {T}
    @assert length(xs) == length(ys) == c.N + c.D + 2
    R = c.ring
    z = gen(R)
    dpi = c.dense_polynomial_interpolator
    # F is the polynomial of degree at max N + D + 1
    # that passes through interpolation points
    F = interpolate!(dpi, xs, ys)
    modulo = prod(z - x for x in xs)
    k = c.N
    r, t, _ = constrainedEEA(F, modulo, k)
    !isunit(gcd(r, t)) && throw("Cauchy interpolation fail.")
    normfactor = trailing_coefficient(t)
    divexact(r, normfactor), divexact(t, normfactor)
end
