
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

function Base.copy(a::adaptiveVanDerHoevenLecerf)
    adaptiveVanDerHoevenLecerf(a.ring, a.N, a.D, a.idx, 
        copy(a.ξa), 
        copy(a.ωi), 
        copy(a.ξij), 
        copy(a.fij),
        copy(a.shift),
        copy(a.univariate_rational_interpolator), 
        copy(a.multivariate_poly_interpolator_num), 
        copy(a.multivariate_poly_interpolator_den)
    )
end

function several_adaptiveVanDerHoevenLecerf(R, N, D, v)
    one_for_all = adaptiveVanDerHoevenLecerf(R, N, D)
    other = [
        [
            adaptiveVanDerHoevenLecerf(R, N, D, shift=one_for_all.shift)
            for i in 1:length(v[j])
        ]
        for j in 1:length(v)
    ]
    other
end

function next_point!(avdhls::Vector{Vector{T}}) where {T}
    a = avdhls[1][1]
    shift = a.shift
    if a.idx == a.N + a.D + 2
        a.idx = 0
    end
    a.idx += 1
    if isone(a.idx)
        ωi = next_point!(a.multivariate_poly_interpolator_num)
        ωi = next_point!(a.multivariate_poly_interpolator_den)
        ξij = distinct_points(base_ring(a.ring), a.N + a.D + 2)
        ωξij = [ωi .* ξ for ξ in ξij]
        # shift in each variable,
        # f(ξ*ωi0 + s0, ξ*ωi1 + s1,ξ*ωi2 + s2,..., ξ*ωin + sn) for ξ in ξij
        ωξsij = [ωξ .+ a.shift for ωξ in ωξij]
        # homogenize
        # ωξsijh = { (ξ*ωi1 + s1)/(ξ*ωi0 + s0), (ξ*ωi2 + s2)/(ξ*ωi0 + s0)... } for ξ in ξij
        ωξsijh = [(ωξs .// ωξs[1])[2:end] for ωξs in ωξsij]
        a.ξa = ωξsijh
        a.ξij = ξij
        a.ωi = ωi
        for i in 1:length(avdhls)
            for j in 1:length(avdhls[i])
                if i == 1 && j == 1
                    continue
                end
                aij = avdhls[i][j]
                @assert shift == aij.shift
                next_point!(aij.multivariate_poly_interpolator_num)
                next_point!(aij.multivariate_poly_interpolator_den)
                aij.ξa = a.ξa
                aij.ξij = a.ξij
                aij.ωi = a.ωi
                aij.idx = a.idx
            end
        end
    end
    a.ξa[a.idx]
end

function next_point!(avdhl::adaptiveVanDerHoevenLecerf)
    next_point!([[avdhl]])
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
        empty!(uri)
        @assert isone(trailing_coefficient(Q))
        inum = avdhl.multivariate_poly_interpolator_num
        iden = avdhl.multivariate_poly_interpolator_den
        cnum = leading_coefficient(P)
        cden = leading_coefficient(Q)
        success_num, num = next!(inum, cnum)
        success_den, den = next!(iden, cden)
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
