
mutable struct SharedStateAVDHL{T,Ring,UnivRing,UnivPoly,Poly}
    ring::Ring
    sjm::SimultaneousJavadiMonagan{T,Ring,UnivRing,UnivPoly,Poly}
    shift::Vector{T}
    ξij::Vector{T}
    N::Int
    D::Int
    idx::Int
end

function SharedStateAVDHL(ring, sjm::SimultaneousJavadiMonagan{T}, N::Integer, D::Integer, shift::Vector{T}) where {T}
    @assert N >= 0 && D >= 0
    SharedStateAVDHL(
        ring, 
        sjm, 
        shift, 
        Vector{T}(undef, N + D + 2), 
        N, D, 
        0
    )
end

mutable struct AdaptiveVanDerHoevenLecerf{T,Ring,UnivRing,UnivPoly,Poly,I1<:AbstractRationalInterpolator} <: AbstractRationalInterpolator
    shared::SharedStateAVDHL{T,Ring,UnivRing,UnivPoly,Poly}
    univariate_rational_interpolator::I1
    fij::Vector{T}
end

function AdaptiveVanDerHoevenLecerf(shared::SharedStateAVDHL{T,Ring,UnivRing}) where {T,Ring,UnivRing}
    N, D = shared.N, shared.D
    univariate_rational_interpolator = Cauchy(univariatize(Nemo.PolyRing, shared.ring), N, D)
    AdaptiveVanDerHoevenLecerf(
        shared,
        univariate_rational_interpolator,
        Vector{T}(undef, N + D + 2)
    )
end

mutable struct SimultaneousAdaptiveVanDerHoevenLecerf{T,Ring,UnivRing,UnivPoly,Poly,I1<:AbstractRationalInterpolator} <: AbstractRationalInterpolator
    # a shared state of all interpolators
    shared::SharedStateAVDHL{T,Ring,UnivRing}
    # a vector of interpolators
    avdhls::Vector{AdaptiveVanDerHoevenLecerf{T, Ring, UnivRing, UnivPoly, Poly, I1}}
    # a vector of results as tuples (is interpolated?, numerator, denominator)
    results::Vector{Tuple{Bool, Poly, Poly}}
end

function SimultaneousAdaptiveVanDerHoevenLecerf(ring, count::Integer, N::Integer, D::Integer)
    @assert count > 0
    @assert N >= 0 && D >= 0
    ringhom = homogenize(ring)
    shift = random_point(ringhom)
    sharedjm = SharedStateJM(ringhom)
    degrees = [isodd(i) ? N : D for i in 1:2*count]
    sjm = SimultaneousJavadiMonagan(sharedjm, degrees)
    shared = SharedStateAVDHL(ring, sjm, N, D, shift)
    avdhls = map(_ -> AdaptiveVanDerHoevenLecerf(shared), 1:count)
    results = map(_ -> (false, zero(ring), zero(ring)), 1:count)
    SimultaneousAdaptiveVanDerHoevenLecerf(shared, avdhls, results)
end

function allready(savdhl::SimultaneousAdaptiveVanDerHoevenLecerf)
    allready(savdhl.shared.sjm) 
end

function getresult(savdhl::SimultaneousAdaptiveVanDerHoevenLecerf, i::Integer)
    savdhl.results[i]
end

function nextpoints!(savdhl::SimultaneousAdaptiveVanDerHoevenLecerf)
    nextpoints!(savdhl.shared)
end

function nextevaluations!(savdhl::SimultaneousAdaptiveVanDerHoevenLecerf{T}, vs) where {T}
    vals = Vector{T}(undef, 2*length(vs))
    for i in 1:length(vs)
        nextevaluations!(savdhl.avdhls[i], i, vals, vs[i])
    end
    nextevaluation!(savdhl.shared.sjm, vals)
    for i in 1:length(vs)
        update!(savdhl.avdhls[i], i, savdhl.results)
    end
    savdhl.results
end

function nextpoints!(shared::SharedStateAVDHL)
    K = base_ring(shared.ring)
    N, D = shared.N, shared.D
    shift = shared.shift
    ωi = nextpoint!(shared.sjm)
    # random points for dense rational interpolation,
    # f(ξ*x0,ξ*x1,ξ*x2,..., ξ*xn) for ξ in ξij
    ξij = distinct_points(K, N + D + 2)
    # "substitute" ωi,
    # f(ξ*ωi0, ξ*ωi1,ξ*ωi2,..., ξ*ωin) for ξ in ξij
    ωξij = [ωi .* ξ for ξ in ξij]
    # shift in each variable,
    # f(ξ*ωi0 + s0, ξ*ωi1 + s1,ξ*ωi2 + s2,..., ξ*ωin + sn) for ξ in ξij
    ωξsij = [ωξ .+ shift for ωξ in ωξij]
    # homogenize
    # ωξsijh = { (ξ*ωi1 + s1)/(ξ*ωi0 + s0), (ξ*ωi2 + s2)/(ξ*ωi0 + s0)... } for ξ in ξij
    ωξsijh = [(ωξs .// ωξs[1])[2:end] for ωξs in ωξsij]
    shared.ξij = ξij
    # @info "nextpoints!" N D ωi
    ωξsijh
end

function nextevaluations!(avdhl::AdaptiveVanDerHoevenLecerf, i, vals, fij)
    # evaluate the blackbox
    # multiply by the correcton factor ω0^(N - D)
    N, D = avdhl.shared.N, avdhl.shared.D
    shift = avdhl.shared.shift
    ω0 = avdhl.shared.sjm.shared.β[1]
    ξij = avdhl.shared.ξij
    fij = map(cξ -> cξ[1] * (ω0*cξ[2] + shift[1])^(N - D), zip(fij, ξij))
    # interpolate the numerator and the denominator densely.
    P, Q = interpolate!(avdhl.univariate_rational_interpolator, ξij, fij)
    @assert isone(trailing_coefficient(Q))
    empty!(avdhl.univariate_rational_interpolator)
    # n and d are the true dense degrees of the numerator and denominator
    n, d = degree(P), degree(Q)
    Pc = leading_coefficient(P)
    Qc = leading_coefficient(Q)
    vals[2i-1] = Pc
    vals[2i] = Qc
    # @info "nextevaluations!" P Q ω0
    nothing
end

function update!(avdhl::AdaptiveVanDerHoevenLecerf, i, results)
    if isready(avdhl.shared.sjm, 2i-1) && isready(avdhl.shared.sjm, 2i)
        _, P = getresult(avdhl.shared.sjm, 2i-1)
        _, Q = getresult(avdhl.shared.sjm, 2i)
        R = avdhl.shared.ring
        xs0 = gens(R)
        xs0 = [one(R), xs0...]
        num = evaluate(P, xs0)
        den = evaluate(Q, xs0)
        normalization_factor = leading_coefficient(den)
        num = map_coefficients(c -> div(c, normalization_factor), num)
        den = map_coefficients(c -> div(c, normalization_factor), den)
        # @warn "ready!" i
        results[i] = true, num, den
    end
    nothing
end

# function nextpoint!(avdhls::Vector{Vector{T}}) where {T}
#     a = avdhls[1][1]
#     shift = a.shift
#     if a.idx == a.N + a.D + 2
#         a.idx = 0
#     end
#     a.idx += 1
#     if isone(a.idx)
#         ωi = next_point!(a.multivariate_poly_interpolator_num)
#         ωi = next_point!(a.multivariate_poly_interpolator_den)
#         ξij = distinct_points(base_ring(a.ring), a.N + a.D + 2)
#         ωξij = [ωi .* ξ for ξ in ξij]
#         # shift in each variable,
#         # f(ξ*ωi0 + s0, ξ*ωi1 + s1,ξ*ωi2 + s2,..., ξ*ωin + sn) for ξ in ξij
#         ωξsij = [ωξ .+ a.shift for ωξ in ωξij]
#         # homogenize
#         # ωξsijh = { (ξ*ωi1 + s1)/(ξ*ωi0 + s0), (ξ*ωi2 + s2)/(ξ*ωi0 + s0)... } for ξ in ξij
#         ωξsijh = [(ωξs .// ωξs[1])[2:end] for ωξs in ωξsij]
#         a.ξa = ωξsijh
#         a.ξij = ξij
#         a.ωi = ωi
#         for i in 1:length(avdhls)
#             for j in 1:length(avdhls[i])
#                 if i == 1 && j == 1
#                     continue
#                 end
#                 aij = avdhls[i][j]
#                 @assert shift == aij.shift
#                 next_point!(aij.multivariate_poly_interpolator_num)
#                 next_point!(aij.multivariate_poly_interpolator_den)
#                 aij.ξa = a.ξa
#                 aij.ξij = a.ξij
#                 aij.ωi = a.ωi
#                 aij.idx = a.idx
#             end
#         end
#     end
#     a.ξa[a.idx]
# end

# function next_point!(avdhl::adaptiveVanDerHoevenLecerf)
#     next_point!([[avdhl]])
# end

# function next!(avdhl::adaptiveVanDerHoevenLecerf, y)
#     R = avdhl.ring
#     N, D = avdhl.N, avdhl.D
#     push!(avdhl.fij, y)
#     if length(avdhl.fij) == N + D + 2
#         # homogenizing component ω0
#         ω0 = avdhl.ωi[1]
#         fij = map(cξ -> cξ[1] * (ω0*cξ[2] + avdhl.shift[1])^(N - D), zip(avdhl.fij, avdhl.ξij))
#         uri = avdhl.univariate_rational_interpolator
#         P, Q = interpolate!(uri, avdhl.ξij, fij)
#         empty!(avdhl.fij)
#         empty!(uri)
#         @assert isone(trailing_coefficient(Q))
#         inum = avdhl.multivariate_poly_interpolator_num
#         iden = avdhl.multivariate_poly_interpolator_den
#         cnum = leading_coefficient(P)
#         cden = leading_coefficient(Q)
#         success_num, num = next!(inum, cnum)
#         success_den, den = next!(iden, cden)
#         if success_num && success_den
#             R = avdhl.ring
#             xs0 = gens(R)
#             xs0 = [one(R), xs0...]
#             num = evaluate(num, xs0)
#             den = evaluate(den, xs0)
#             normalization_factor = trailing_coefficient(den)
#             num = map_coefficients(c -> div(c, normalization_factor), num)
#             den = map_coefficients(c -> div(c, normalization_factor), den)
#             return true, (num, den)
#         else
#             return false, (num, den)
#         end
#     end
#     return false, (one(R), one(R))
# end
