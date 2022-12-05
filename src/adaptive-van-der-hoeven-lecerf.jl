
mutable struct adaptiveVanDerHoevenLecerf{Ring, T, I1<:AbstractRationalInterpolator, I2<:AbstractPolynomialInterpolator} <: AbstractRationalInterpolator
    ring::Ring
    N::Int
    D::Int
    idx::Int
    ξa::Vector{Vector{T}}
    ωi::Vector{T}
    ξij::Vector{T}
    fij::Vector{T}
    shift::Vector{T}
    univariate_rational_interpolator::I1
    multivariate_poly_interpolator_num::I2
    multivariate_poly_interpolator_den::I2
end

function adaptiveVanDerHoevenLecerf(ring, N::Integer, D::Integer; 
        shift=random_point(homogenize(ring)),
        univariate_rational_interpolator=Cauchy(univariatize(Nemo.PolyRing, ring), N, D),
        multivariate_poly_interpolator=BenOrTiwari(homogenize(ring)))
    @assert N >= 0 && D >= 0
    adaptiveVanDerHoevenLecerf(ring, N, D, 0, 
        Vector{Vector{elem_type(base_ring(ring))}}(undef, 0), 
        Vector{elem_type(base_ring(ring))}(undef, 0), 
        Vector{elem_type(base_ring(ring))}(undef, 0), 
        Vector{elem_type(base_ring(ring))}(undef, 0),
        shift,
        univariate_rational_interpolator, 
        multivariate_poly_interpolator, 
        copy(multivariate_poly_interpolator)
    )
end

function next_point!(avdhl::adaptiveVanDerHoevenLecerf)
    if iszero(avdhl.idx)
        ωi = next_point!(avdhl.multivariate_poly_interpolator_num)
        # random points for dense rational interpolation,
        # f(ξ*x0,ξ*x1,ξ*x2,..., ξ*xn) for ξ in ξij
        ξij = distinct_points(base_ring(avdhl.ring), avdhl.N + avdhl.D + 2)
        # "substitute" ωi,
        # f(ξ*ωi0, ξ*ωi1,ξ*ωi2,..., ξ*ωin) for ξ in ξij
        ωξij = [ωi .* ξ for ξ in ξij]
        # shift in each variable,
        # f(ξ*ωi0 + s0, ξ*ωi1 + s1,ξ*ωi2 + s2,..., ξ*ωin + sn) for ξ in ξij
        ωξsij = [ωξ .+ avdhl.shift for ωξ in ωξij]
        # homogenize
        # ωξsijh = { (ξ*ωi1 + s1)/(ξ*ωi0 + s0), (ξ*ωi2 + s2)/(ξ*ωi0 + s0)... } for ξ in ξij
        ωξsijh = [(ωξs .// ωξs[1])[2:end] for ωξs in ωξsij]
        avdhl.ξa = ωξsijh
        avdhl.ξij = ξij
        avdhl.ωi = ωi
    end
    avdhl.idx += 1
    ans = avdhl.ξa[avdhl.idx]
    if avdhl.idx == avdhl.N + avdhl.D + 2
        avdhl.idx = 0
    end
    ans
end

function next!(avdhl::adaptiveVanDerHoevenLecerf, y)
    R = avdhl.ring
    N, D = avdhl.N, avdhl.D
    push!(avdhl.fij, y)
    if length(avdhl.fij) == N + D + 2
        # homogenizing component ω0
        ω0 = avdhl.ωi[1]
        fij = map(cξ -> cξ[1] * (ω0*cξ[2] + avdhl.shift[1])^(N - D), zip(avdhl.fij, avdhl.ξij))
        uri = avdhl.univariate_rational_interpolator
        P, Q = interpolate!(uri, avdhl.ξij, fij)
        empty!(avdhl.fij)
        @assert isone(trailing_coefficient(Q))
        empty!(uri)
        inum = avdhl.multivariate_poly_interpolator_num
        iden = avdhl.multivariate_poly_interpolator_den
        cnum = leading_coefficient(P)
        cden = leading_coefficient(Q)
        success_num, num = next!(inum, next_point!(inum), cnum)
        success_den, den = next!(iden, next_point!(iden), cden)
        if success_num && success_den
            R = avdhl.ring
            xs0 = gens(R)
            xs0 = [one(R), xs0...]
            num = evaluate(num, xs0)
            den = evaluate(den, xs0)
            normalization_factor = trailing_coefficient(den)
            num = map_coefficients(c -> div(c, normalization_factor), num)
            den = map_coefficients(c -> div(c, normalization_factor), den)
            return true, (num, den)
        else
            return false, (num, den)
        end
    end
    return false, (one(R), one(R))
end
