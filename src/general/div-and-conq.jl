
# Product of all (z - x) for x in xs[i..j] 
function _producttree(z, xs, i, j)
    i == j && return z - xs[i]
    m = div(i + j, 2)
    _producttree(z, xs, i, m) * _producttree(z, xs, m + 1, j)
end

# Gcd of all f in f[i..j] 
function _gcdtree(f, i, j)
    i == j && return f[i]
    m = div(i + j, 2)
    gcd(_gcdtree(f, i, m), _gcdtree(f, m + 1, j))
end

function gcdtree(f)
    _gcdtree(f, 1, length(f))
end

# computes the product of all (z - x) for x in xs
# using the tree product algorithm.
# O(M(n)), n = length(xs)
function producttree(z, xs)
    _producttree(z, xs, 1, length(xs))
end

treeroot(tree) = tree[end][1]
treedepth(tree) = length(tree)
treebase(tree) = length(first(tree))

# Build the tree of products of (z - x) for x in xs and returns it.
# Leaves of the tree are (z - x) and the root is the product,
# O(M(n)), n = length(xs)
function buildproducttree(z::T, xs) where {T}
    n = length(xs)
    @assert n > 0
    npow = nextpow(2, n)
    k = round(Int, log(2, npow)) + 1
    tree = Vector{Vector{T}}(undef, k)
    tree[1] = Vector{T}(undef, npow)
    @inbounds for i in 1:n
        tree[1][i] = z - xs[i]
    end
    @inbounds for i in n+1:npow
        tree[1][i] = one(z)
    end
    @inbounds for i in 2:k
        nel = 2^(k - i)
        tree[i] = Vector{T}(undef, nel)
        for j in 1:nel
            tree[i][j] = tree[i-1][2*j-1] * tree[i-1][2*j]
        end
    end
    tree
end

function _remindertree!(f, rtree, ptree, depth, idx)
    iszero(f) && return rtree
    if depth == 0
        rtree[idx] = coeff(mod(f, ptree[1][idx]), 0)
        return rtree
    end
    l, r = 2 * idx - 1, 2 * idx
    @inbounds r0 = mod(f, ptree[depth][l])
    @inbounds r1 = mod(f, ptree[depth][r])
    _remindertree!(r0, rtree, ptree, depth - 1, l)
    _remindertree!(r1, rtree, ptree, depth - 1, r)
    rtree
end

# Given a polynomial f and a product tree ptree
# built from (z - x) for x in xs,
# returns the vector of reminders (f mod (z - x)) for x in xs.
# O(M(n)logn), if n = degree(f) = O(length(xs))
function remindertree(f, ptree)
    K = base_ring(parent(f))
    rtree = zeros(K, treebase(ptree))
    _remindertree!(f, rtree, ptree, treedepth(ptree) - 1, 1)
end

# Solves the system
# | 1       1     ...     1    | | x1 |   | a1 |
# | v1      v2    ...     vT   | | x2 |   | a2 |
# | .       .     ...     .    | | .  | = | .  |
# | vn^k    vn^k  ...     vn^k | | xn |   | an |
# where k is n-1
# O(M(n)logn)
# Notations are taken from 
#   "Improved Sparse Multivariate Polynomial Interpolation Algorithms",
#   Erich Kaltofen and Lakshman Yagati
function solve_transposed_vandermonde(Rz, vi, ai)
    # @assert length(vi) == length(ai)
    t, T = length(vi), length(ai)
    z = gen(Rz)
    # Dz = a1*z^n + a2*z^n-1 + ... + an*z
    # O(1)
    Dz = Nemo.shift_left(Rz(reverse(ai)), 1)
    # Bz = (z - v1)(z - v2)...(z - vn)
    # O(M(n))
    ptree = buildproducttree(z, vi)
    Bz = treeroot(ptree)
    ∂B = derivative(Bz)
    # αi = ∏(vi - vj) for all j ≠ i
    # O(M(n)logn)
    αi = remindertree(∂B, ptree)
    # O(M(n))
    Qn = Nemo.shift_right(Bz * Dz, t + 1)
    # Qnvi = Qn evaluated at vi for all i
    # O(M(n)logn)
    Qnvi = remindertree(Qn, ptree)
    # O(n)
    xi = map(i -> Qnvi[i] * inv(αi[i]), 1:t)
    xi
end

function _solve_toeplitz(R, n, x, y, b::Vector{T}) where {T}
    t = gen(R)
    m = t^(n + 1)
    X, Y = x, y
    B = R(b)
    P = mod(reverse(X, n + 2) * reverse(B, n + 1), m)
    Q = mod(Y * reverse(B, n + 1), m)
    Z = mod(X * reverse(Q, n + 1) - reverse(Y, n + 2) * reverse(P, n + 1), m)
    Z = div(Z, trailing_coefficient(X))
    Z
end

function first_row_and_col_of_toeplitz_inverse(R, n::Int, a::AbstractVector)
    K, t = base_ring(R), gen(R)
    U0 = t^(2n + 1)
    U1 = R(a)
    Uj, Wj, Vj = fastconstrainedEEA(U0, U1, n)
    # the system is singular
    degree(Uj) < n && return false, a, a
    @assert !iszero(trailing_coefficient(Vj))
    x = inv(leading_coefficient(Uj)) * Vj
    U1 = reverse(U1)
    Uj, Wj, Vj = fastconstrainedEEA(U0, U1, n)
    y = inv(leading_coefficient(Uj)) * Vj
    return true, x, y
    # NOTE: the original algorithm includes MORE steps.
    # We limit outselves to the easy case and implement less steps
end

#=
    Solves the Toeplitz system
    Ax = b
=#
function solve_toeplitz(R, a::Vector{T}, b::Vector{T}) where {T}
    @assert !isempty(a) && !isempty(b)
    n = length(b) - 1
    @assert length(a) == 2n + 1
    solvable, x, y = first_row_and_col_of_toeplitz_inverse(R, n, a)
    @assert solvable
    _solve_toeplitz(R, n, x, y, b)
end

function _lagrangetree(z, ys, ptree, depth, idx)
    R = parent(z)
    if depth == 0
        return R(ys[idx])
    end
    l, r = 2 * idx - 1, 2 * idx
    r0 = _lagrangetree(z, ys, ptree, depth - 1, l)
    r1 = _lagrangetree(z, ys, ptree, depth - 1, r)
    r0 * ptree[depth][r] + r1 * ptree[depth][l]
end

function lagrangetree(z, ys, ptree)
    _lagrangetree(z, ys, ptree, treedepth(ptree) - 1, 1)
end

# Returns a unique univariate polynomial f in the ring R,
# such that f(x) = y for all (x, y) in xs, ys;
# O(M(n)logn), where n = O(length(xs))
function fastpolyinterpolate(R, xs, ys)
    @assert length(xs) == length(ys)
    z = gen(R)
    # O(M(n))
    ptree = buildproducttree(z, xs)
    m = treeroot(ptree)
    dm = derivative(m)
    # O(M(n)logn)
    si = remindertree(dm, ptree)
    ysi = zeros(base_ring(R), nextpow(2, length(ys)))
    for i in 1:length(ys)
        ysi[i] = ys[i] * inv(si[i])
    end
    # O(M(n)logn)
    lagrangetree(z, ysi, ptree)
end


# Adapted from
#   "Solving Systems of Non-Linear Polynomial Equations Faster"
#   by Canny, Kaltofen, and Yagati, 1989
# and
#   "Fast Solution of Toeplitz Systems of Equations and 
#   Computation of Pade Approximants" 
#   by Brent, Gustavson, Yun, 1980
#
# Given a multivariate polynomial
#   f(x1..xn) = a1m1 + a2m2 + ... + atmt, 
# an evaluation point
#   ω = [ω1^i,...,ωn]
# Computes and returns b = (b1,..,bt) = (f(ω^0),..,f(ω^{t-1}))
# O(M(t)log t) arithmetic operations in K.
#
# !!! This is internal function.
#     It assumes that the number of evaluation points T
#     is equal to or greater than the number of monomials t.
function _fast_multivariate_evaluate(R, f, ω, T)
    # vi = mi(ω)
    #     |1 v1 ... v1^{T-1}|
    # V = |1 v2 ... v2^{T-1}|
    #     |     ...
    #     |1 vt ... vt^{T-1}|
    # a = (a1...at), b = (b0...bT-1)
    # Then
    #   b = V^T a  (1)
    # Let a' = V \ a (solve this in M(t)log t)
    # Then (1) becomes 
    #   b = (V^T V) a'
    # The product of Hankel matrix (V^T V) and a'
    # can be computed in O(M(t)), and the entries of
    # V^T V can be computed in O(M(T)log t)
    n = nvars(R)
    K = base_ring(R)
    Runiv, z = PolynomialRing(K, "z")
    @assert n > 1
    a = collect(coefficients(f))
    t = length(a)
    @assert T >= T

    # O(M(T)log t)
    vi = map(i -> evaluate(monomial(f, i), ω), 1:t)
    a_prime = fastpolyinterpolate(Runiv, vi, a)

    # sigmas are elementary symmetric functions of vi:
    # sigma1 = v1 +..+ vn, 
    # sigma2 = v1v2 + v1v2 +..,
    # sigmat = v1v2..vt 
    sigmai = reverse(collect(coefficients(producttree(z, vi)))[1:t])
    Ai = vcat([zero(K) for _ in 1:2T-1], [one(K)], sigmai, [zero(K) for _ in 1:(2T-t-1)])
    wi = vcat(map(i -> (-i) * sigmai[i], 1:t), [zero(K) for _ in 1:2T-t])

    # sj are the power sums of vi:
    # s1 = v1 +..+ vt,
    # s2 = v1^2 +..+ vt^2,
    # st = v1^t +..+ vt^t
    sj = solve_toeplitz(Runiv, Ai, wi)
    sj = shift_left(sj, 1) + t

    g = Runiv(reverse(a_prime, t))
    res = sj * g
    res = mod(shift_right(res, t - 1), z^(T))

    collect(coefficients(res))
end

# Dispatch between the two cases, t = # term in f:
# t <= T and t > T
function fast_multivariate_evaluate(R, f, ω, T::Integer)
    t = length(f)
    K = base_ring(R)
    @assert T > 0
    t <= T && return _fast_multivariate_evaluate(R, f, ω, T)
    # here t > T
    q = div(t, T)
    evals = map(_ -> zero(K), 1:T)
    tms = collect(terms(f))
    for i in 1:q
        fi = sum(view(tms, (1+T*(i-1)):T*i))
        evi = _fast_multivariate_evaluate(R, fi, ω, T)
        @assert length(evi) == T
        for j in 1:T
            evals[j] += evi[j]
        end
    end
    rm = t - q * T
    if !iszero(rm)
        fend = sum(tms[T*q+1:end])
        evend = _fast_multivariate_evaluate(R, fend, ω, T)
        for j in 1:T
            evals[j] += evend[j]
        end
    end
    evals
end

# Given a vector of multivariate polynomials {F_j}, j = 1..k and an evaluation
# point, returns the array of evaluations: {F_j(p^i)}, j = 1..k, i = 0..T-1
function simultaneous_multivariate_evaluate(R, polys, point, T)
    field = base_ring(R)
    polys_eval = Vector{Vector{eltype(point)}}(undef, T)
    for i in 1:length(polys_eval)
        polys_eval[i] = map(_ -> zero(field), 1:length(polys))
    end
    for i in 1:length(polys)
        iszero(polys[i]) && continue
        coeffs_eval = fast_multivariate_evaluate(R, polys[i], point, T)
        for j in 1:length(coeffs_eval)
            polys_eval[j][i] = coeffs_eval[j]
        end
    end
    polys_eval
end
