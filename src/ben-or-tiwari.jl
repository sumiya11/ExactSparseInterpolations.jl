#=
    Ben-Or and Tiwari algorithm.

    Univariate/Multivariate, Sparse, Deterministic, 
    Blackbox
    
    2*T evaluations, where T is the upper bound on the number of terms.

    The algorithm comes in two flavors:
    Deterministic explicit, and MC Probabilistic adaptive

    Explicit version:
        2*T queries to blackbox, T is the upper bound on t
    Adaptive version:
        2*t queries to blackbox, t is the number of terms
=#

#=
    !!!
    Any nonadaptive algorithm 
    for sparse polynomial interpolation 
    must perform 2t evaluations in the worst case.
    In that sence, Ben-Or-Tiwari is optimal.
=#

mutable struct BenOrTiwari{T1,T2,T3,T4} <: AbstractInterpolator
    ring::T1
    t::Int
    T::Int
    ω::Vector{T2}
    vi::Vector{T2}
    i::Int
    bm::BerlekampMassey{T2, T3, T4}
end

function BenOrTiwari(ring; 
                t::Integer=-1, T::Integer=-1, 
                ω=geometric_point(ring),
                bm=BerlekampMassey(true_univariatize(ring)))
    BenOrTiwari(ring, t, T, ω, empty(ω), 1, bm)
end

function integer_to_monom(m, u0)
    d = Primes.factor(Dict, m)
    monom = zeros(Int, length(u0))
    for (i, p) in enumerate(u0)
        p = BigInt(numerator(p))
        if haskey(d, p)
            monom[i] = d[p]
        end
    end
    monom
end

function next!(bot::BenOrTiwari{T1,T2,T3,T4}, x, y) where {T1,T2,T3,T4}
    R = parent(bot)
    K = base_ring(R)
    ω = bot.ω
    i = bot.i
    ωi = ω .^ i
    @assert ωi == x
    vi = y
    bm = bot.bm
    bot.i += 1
    push!(bot.vi, y)
    success, L = next!(bm, vi)
    # @info "" success, L
    if success
        t = degree(reverse(L))
        mi = Nemo.roots(reverse(L))
        # @info "" mi
        monoms = [
            integer_to_monom(BigInt(numerator(m)), ω)
            for m in mi
        ]
        # @info "" monoms
        Vspace = Nemo.MatrixSpace(K, t, t)
        A = Vspace(reshape([
            mi[i] ^ (j)
            for i in 1:t
                for j in 1:t
        ], (t, t)))
        vi = bot.vi[1:t]
        ai = Nemo.inv(A) * vi
        return (success, R(ai, monoms))
    else
        return (success, one(R)) 
    end
end

function Nemo.interpolate(bot::BenOrTiwari{T1}, blackbox) where {T1}
    R = parent(bot)
    K = base_ring(R)
    ω = bot.ω
    i = bot.i
    success = false
    f = one(R)
    while !success
        ωi = ω .^ i
        vi = blackbox(ωi)
        success, f = next!(bot, ωi, vi)
        i += 1
    end
    f
end

# function Nemo.interpolate(bot::BenOrTiwari{T1}, blackbox) where {T1}
#     R = parent(bot)
#     K = base_ring(R)
#     u0 = geometric_point(R)
#     t = bot.t
#     ui = [u0 .^ i for i in 0:2t-1]
#     vi = map(blackbox, ui)
#     Vspace = Nemo.MatrixSpace(K, t, t)
#     V = Vspace(reshape([
#         vi[i + j - 1]
#         for i in 1:t
#             for j in 1:t
#     ], (t, t)))
#     @info "" ui vi V
#     # k = Nemo.rank(V)
#     # s = [vi[i + k] for i in 1:t]
#     # Vbar = Vspace(reshape([
#     #     vi[j - i + 1]
#     #     for i in 1:t
#     #         for j in 1:t
#     # ], (t, t)))
#     # lambda = Nemo.solve(Vbar, s)
#     k = t
#     s = [vi[i+k] for i in 1:t]
#     lambda = Nemo.inv(V) * s
#     @info "" s lambda
#     Ru, z = PolynomialRing(K, "z")
#     L = z^k - sum([lambda[i+1]*z^i for i in 0:k-1])
#     @info "" L
#     mi = Nemo.roots(L)
#     @info "" mi
#     monoms = [
#         integer_to_monom(BigInt(numerator(m)), u0)
#         for m in mi
#     ]
#     @info "" monoms
#     A = Vspace(reshape([
#         mi[i] ^ (j-1)
#         for i in 1:t
#             for j in 1:t
#     ], (t, t)))
#     vi = vi[1:t]
#     ai = Nemo.inv(A) * vi
#     @info "" ai
#     R(ai, monoms)
# end
