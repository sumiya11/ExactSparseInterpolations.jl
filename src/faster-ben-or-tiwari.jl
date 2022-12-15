
# Interpolate polynomials in O(M(n)logn),
# where M(n) is the complexity of polynomial multiplication
mutable struct FasterBenOrTiwari{Ring} <: AbstractPolynomialInterpolator
    # multivariate or univariate polynomial ring
    # (currently only univariate ones are supported)
    ring::Ring
    # the number of terms in the interpolant
    T::Int
end

# Given (polynomials) f and g (|f| <= |g|),
# computes and returns a single row from the EEA algorithm (r, t, s), 
# such that r = t*f + s*g, |r| < k, where |r| is maximal possible.
# O(M(T)logT), where T = max(degree(f), degree(g))
function Padé(f, g, k)
    r, t, s = fastconstrainedEEA(g, f, k)
    r, s, t
end

# Solve the system
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
    @assert length(vi) == length(ai)
    n = length(vi)
    z = gen(Rz)
    # Dz = a1^z^n + a2^z^n-1 + ... + an^z
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
    Qn = Nemo.shift_right(Bz*Dz, n + 1)
    # Qnvi = Qn evaluated at vi for all i
    # O(M(n)logn)
    Qnvi = remindertree(Qn, ptree)
    # O(n)
    xi = map(i -> Qnvi[i] * inv(αi[i]), 1:n)
    xi
end

# Returns an appropriate starting point 
# for the geometric sequence ω^0, ω^1, ω^2, ...
function startingpoint(fbot::FasterBenOrTiwari)
    R = fbot.ring
    K = base_ring(R)
    # if univariate
    if nvars(R) == 1
        if degree(K) == 1
            # if the field is actually Z/Zp,
            # return an element of order p-1
            [randomgenerator(K)]
        else
            # otherwise, return a primitive element
            [gen(K)]
        end
    else
        throw("Not implemented yet, sorry!")
    end
end

function interpolate!(fbot::FasterBenOrTiwari, blackbox)
    Rx = fbot.ring
    K = base_ring(Rx)
    T = fbot.T
    # the base of the geometric sequence ω^1, ω^2, ...
    ω = startingpoint(fbot)
    # generate the sequence
    # O(TlogT)
    ωs = map(i -> ω.^i, 0:2T-1)
    # evaluate the blackbox function at the sequence
    # O(LT), where L is the cost of evaluating the blackbox function
    ys = map(blackbox, ωs)
    # construct the polynomial ys[1]z^1 + ys[2]z^2 + ... + ys[2T]z^(2T)
    # O(1)
    Rz, z = K["z"]
    sequence = Rz(ys)
    # find A/B such that A/B = sequence mod z^(2T) and degree(A) < T
    # O(M(T)logT)
    _, B, _ = Padé(sequence, z^(2T), T-1)
    @assert degree(B) == T
    # assuming this is O(T logT^k logq^m) for some k and m, 
    # where q is the order of the base field
    mi = map(inv, Nemo.roots(B))
    @assert length(mi) == T
    # find the monomials of the interpolant,
    # O(TlogTlogq), where q is the order of the base field
    # (note that this cost covers the case where K is not a prime field)
    # (assuming ord is smooth)
    monoms = map(m -> discrete_log(ω, m, PrecomputedField(K)), mi)
    # find the coefficients of the interpolant
    # by solving a T×T Vandermonde system
    # O(M(T)logT)
    coeffs = solve_transposed_vandermonde(Rz, view(mi, 1:T), view(ys, 1:T))
    Rx(coeffs, map(e -> map(Int, e), monoms))
end
