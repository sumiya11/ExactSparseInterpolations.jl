#=
    Ben-Or and Tiwari algorithm.

    Univariate/Multivariate, Sparse, Deterministic, Blackbox
    
    2*T evaluations, where T is the upper bound on the number of terms t.

    The algorithm comes in two flavors:
    Deterministic, and MC Probabilistic Adaptive

    Deterministic version:
        2*T queries to the blackbox
    Probabilistic version:
        2*t queries to the blackbox
=#

#=
    Any non-adaptive algorithm 
    for sparse polynomial interpolation 
    must perform 2t evaluations in the worst case.
    In that sence, Ben-Or-Tiwari is optimal.
=#

mutable struct BenOrTiwari{Ring,K,UnivRing,UnivPoly} <: AbstractPolynomialInterpolator
    # multivariate polynomial ring
    ring::Ring
    # exact number of terms in the interpolant
    t::Int
    # upper bound on the number of terms in the interpolant
    T::Int
    # initial evaluation point
    ω0::Vector{K}
    # current evaluation point
    ω::Vector{K}
    # vector of blackbox values at ωi (initially empty)
    vi::Vector{K}
    i::Int
    bm::BerlekampMassey{UnivRing, K, UnivPoly}
end

function BenOrTiwari(ring::Ring; 
                t::Union{Integer, Nothing}=nothing, 
                T::Union{Integer, Nothing}=nothing,
                bm=BerlekampMassey(univariatize(AbstractAlgebra.PolyRing, ring))) where {Ring}
    # check that input parameters are okay
    @assert isnothing(t) || isnothing(T) || t <= T
    isnothing(T) && (T = typemax(Int))
    isnothing(t) && (t = min(typemax(Int), T))
    ω0 = zeros(base_ring(ring), 0)
    ω = zeros(base_ring(ring), 0)
    vi = zeros(base_ring(ring), 0)
    BenOrTiwari(ring, t, T, ω0, ω, vi, 0, bm)
end

function Base.copy(bot::BenOrTiwari{Ring,K,UnivRing,UnivPoly}) where {Ring,K,UnivRing,UnivPoly}
    BenOrTiwari{Ring,K,UnivRing,UnivPoly}(bot.ring, bot.t, bot.T, copy(bot.ω0), copy(bot.ω), copy(bot.vi), bot.i, copy(bot.bm))
end

# Returns the next point for blackbox evaluation
function next_point!(bot::BenOrTiwari; increment=false)
    nv = nvars(bot.ring)
    K = base_ring(bot.ring)
    # if this is the first point in the sequence
    if iszero(bot.i)
        # ω is [2, 3, 5, ..., prime_nv]
        ps = Primes.nextprimes(2, nv)
        bot.ω0 = map(K, ps)
        bot.ω = bot.ω0
    else
        # ω is [2^i, 3^i, 5^i, ..., prime_nv^i]
        bot.ω = map(x -> x^(bot.i+1), bot.ω0)
    end
    increment && (bot.i += 1)
    bot.ω
end

# Given an integer m and a vector of primes p = [p1,p2,..pk]
# returns a vector [i1,i2,..,ik], s.t. m == p1^i1*p2^i2*...*pk^ik
discrete_log(m, p::AbstractVector) = discrete_log(BigInt(m), p)

function discrete_log(m::BigInt, p::AbstractVector)
    d = Primes.factor(Dict, m)
    monom = zeros(Int, length(p))
    for (i, p) in enumerate(p)
        p = BigInt(numerator(p))
        if haskey(d, p)
            monom[i] = d[p]
        end
    end
    monom
end

# 
function next!(bot::BenOrTiwari{Ring,K,UnivRing,UnivPoly}, v::K) where {Ring,K,UnivRing,UnivPoly}
    R = bot.ring
    bot.i += 1
    push!(bot.vi, v)
    # make one iteration of Berlekamp-Massey for input point v
    success, L = next!(bot.bm, v)
    # if the generator of the sequence was Not found
    !success && return (false, one(R))
    # if L is the generator of the sequence
    Lrev = reverse(L)
    t = degree(Lrev)
    mi = Nemo.roots(Lrev)
    # find the monomials of the interpolant:
    monoms = [
        discrete_log(numerator(m), bot.ω0)
        for m in mi
    ]
    # find the coefficients of the interpolant:
    # matrix A is a t×t Vandermonde matrix,
    # solve A*ai = vi for ai
    Vspace = Nemo.MatrixSpace(base_ring(R), t, t)
    A = Vspace(reshape([
        mi[i] ^ (j)
        for i in 1:t
            for j in 1:t
    ], (t, t)))
    vi = bot.vi[1:t]
    ai = Nemo.inv(A) * vi
    true, R(ai, monoms)
end

function interpolate!(bot::BenOrTiwari, blackbox)
    success = false
    f = one(bot.ring)
    while !success
        x = next_point!(bot)
        y = blackbox(x)
        success, f = next!(bot, y)
    end
    f
end
