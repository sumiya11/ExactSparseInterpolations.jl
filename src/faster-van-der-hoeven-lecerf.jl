
# Measure running time of each step!
const _runtime_vdhl_cauchy = Float64[]
const _runtime_vdhl_pointsmanip = Float64[]
const _runtime_vdhl_evaluate = Float64[]
const _runtime_vdhl_benortiwari = Float64[]

function _runtime_vdhl_dump()
    k = div(length(_runtime_vdhl_benortiwari), 2)
    cauchy = sum(_runtime_vdhl_cauchy) / k
    pointsmanip = sum(_runtime_vdhl_pointsmanip) / k
    evaluate = sum(_runtime_vdhl_evaluate) / k
    benortiwari = sum(_runtime_vdhl_benortiwari) / k
    ans = (
        cauchy=cauchy,
        pointsmanip=pointsmanip,
        evaluate=evaluate,
        benortiwari=benortiwari,
    )
    empty!(_runtime_vdhl_cauchy)
    empty!(_runtime_vdhl_pointsmanip)
    empty!(_runtime_vdhl_evaluate)
    empty!(_runtime_vdhl_benortiwari)
    ans
end

# Faster implementation of van-der-Hoeven & Lecers rational interpolation.
#
# Needs to know the total degrees, partial degrees, and the number of terms
# of numerator and denominator of the interpolant.
mutable struct FasterVanDerHoevenLecerf{Ring,UnivRing}
    # multivariate polynomial ring
    ring::Ring
    # the total degrees of the numerator/denominator
    Nd::Int
    Dd::Int
    # vectors of partial degrees in the numerator/denominator
    Nds::Vector{Int}
    Dds::Vector{Int}
    # the number of terms in the numerator/denominator
    Nt::Int
    Dt::Int
    # dense interpolator for univariate functions
    cauchy::FasterCauchy{UnivRing}
    # polynomial interpolators for the numerator/denominator
    Ni::FasterMultivariateBenOrTiwari{Ring}
    Di::FasterMultivariateBenOrTiwari{Ring}

    # Construct an FasterVanDerHoevenLecerf object from the
    # collected info of the interpolant' degrees/terms 
    function FasterVanDerHoevenLecerf(
        ring::Ring,
        boundsinfo::NamedTuple
        ) where {Ring}
        FasterVanDerHoevenLecerf(
            ring,
            boundsinfo.numtotaldeg, boundsinfo.dentotaldeg,
            boundsinfo.numpartialdegs, boundsinfo.denpartialdegs,
            boundsinfo.numnterms, boundsinfo.dennterms
        )
    end

    function FasterVanDerHoevenLecerf(
        ring::Ring,
        Nd::Int, Dd::Int,
        Nds::Vector{<:Integer}, Dds::Vector{<:Integer},
        Nt::Int, Dt::Int
        ) where {Ring}
        @assert Nd >= 0 && Dd >= 0
        @assert Nt >= 0 && Dt >= 0
        @assert all(>=(0), Nds) && all(>=(0), Dds)
        n = nvars(ring)
        @assert length(Dds) == length(Nds) == n
        K = base_ring(ring)
        Runiv, _ = Nemo.PolynomialRing(K, "u")
        cauchy = FasterCauchy(Runiv, Nd, Dd)
        ringhom = homogenize(ring)
        Nds = [Nd, Nds...]
        Dds = [Dd, Dds...]
        # (!) we take the maximum for each of the partial degrees,
        # in order to be able use the same points for numerator/denominator
        NDds = map(maximum, zip(Nds, Dds))
        Ni = FasterMultivariateBenOrTiwari(ringhom, Nt, NDds)
        Di = FasterMultivariateBenOrTiwari(ringhom, Dt, NDds)
        new{Ring,typeof(Runiv)}(
            ring, 
            Nd, Dd, 
            Nds, Dds, 
            Nt, Dt, 
            cauchy, 
            Ni, Di
        )
    end
end

# Interpolate the blackbox!
# in T*M(D)log(D) + O(T*n*D*log(q)) + T*D*L + M(T)log(T)
function interpolate!(vdhl::FasterVanDerHoevenLecerf, blackbox)
    # Polynomial ring K[x1,x2...xn]
    R = vdhl.ring
    K = base_ring(R)
    Nd, Dd = vdhl.Nd, vdhl.Dd
    Nt, Dt = vdhl.Nt, vdhl.Dt
    Ni, Di = vdhl.Ni, vdhl.Di
    cauchy = vdhl.cauchy
    T = max(Nt, Dt)
    # Polynomial ring K[x0,x1,x2...xn]
    Rhom = Ni.ring
    # We will substitute points, such that
    # point[j]*dilation[j] + shift[j]
    shift = random_point(Rhom)
    dilation = random_point(Rhom)
    # The starting point in the geometric sequence...
    ω = startingpoint(Ni)
    # ... and the sequence itself (the first degree is 0)
    ωs = map(i -> ω .^ i, 0:2T-1)
    Nys = Vector{elem_type(K)}(undef, 2T)
    Dys = Vector{elem_type(K)}(undef, 2T)
    ωξij = Vector{Vector{elem_type(K)}}(undef, Nd + Dd + 2)
    ωξij0 = Vector{elem_type(K)}(undef, Nd + Dd + 2)
    ξij = distinct_points(K, Nd + Dd + 2)
    # This cycle below is 
    # T*((D + 2)*4n*log(q) + D*L + 2*D*log(q) + M(D)log(D)),
    # which is T*M(D)log(D) + O(T*n*D*log(q)) + T*D*L
    for i in 0:2T-1
        ωi = ωs[i + 1]
        ω0 = ωi[1]
        # stats = @timed ...
        # push!(storage, stats.time)
        ξij[1] = random_point(K)
        # The cycle below is (D + 2)*4n*log(q)
        stats = @timed @inbounds for j in 1:Nd + Dd + 2
            !isassigned(ωξij, j) && (ωξij[j] = zeros(K, length(ω) - 1))
            for nj in 2:length(ω)
                ωξij[j][nj - 1] = ωi[nj]
            end
            ξ = ξij[j]
            ωξij0[j] = ω0*ξ*dilation[1] + shift[1]
            for nj in 2:length(ω)
                ωξij[j][nj-1] = ωξij[j][nj-1]*ξ*dilation[nj] + shift[nj]
                ωξij[j][nj-1] = ωξij[j][nj-1] // ωξij0[j]
            end
        end
        push!(_runtime_vdhl_pointsmanip, stats.time)
        # "substitute" ωi,
        # f(ξ*ωi0, ξ*ωi1,ξ*ωi2,..., ξ*ωin) for ξ in ξij
        # ωξij = [ωi .* ξ for ξ in ξij]
        # # shift in each variable,
        # # f(ξ*ωi0 + s0, ξ*ωi1 + s1,ξ*ωi2 + s2,..., ξ*ωin + sn) for ξ in ξij
        # # dilate in each variable,
        # # f(ξ*ωi0*d0 + s0, ξ*ωi1*d0 + s1, ξ*ωi2*d2 + s2,..., ξ*ωin*dn + sn) for ξ in ξij
        # ωξsij = [ωξ .* dilation .+ shift for ωξ in ωξij]
        # # homogenize,
        # # ωξsijh = { (ξ*ωi1*d1 + s1)/(ξ*ωi0*d0 + s0), (ξ*ωi2*d2 + s2)/(ξ*ωi0*d0 + s0)... } for ξ in ξij
        # ωξsijh = [(ωξs .// ωξs[1])[2:end] for ωξs in ωξsij]
        # evaluate the blackbox
        # fij = map(blackbox, ωξsijh)
        # O(D*L)
        stats = @timed fij = map(blackbox, ωξij)
        push!(_runtime_vdhl_evaluate, stats.time)
        # multiply by the correction factor ω0^(N - D)
        # 2*D*log(q)
        stats = @timed fij = map(cξ -> cξ[1] * cξ[2]^(Nd - Dd), zip(fij, ωξij0))
        push!(_runtime_vdhl_pointsmanip, stats.time)
        # interpolate the numerator and the denominator densely.
        # M(D)logD
        stats = @timed N, D = interpolate!(cauchy, ξij, fij)
        push!(_runtime_vdhl_cauchy, stats.time)
        @assert isone(trailing_coefficient(D))
        Nys[i + 1] = leading_coefficient(N)
        Dys[i + 1] = leading_coefficient(D)
    end
    # Interpolate the leading coefficients in the numerator and denominator
    # M(T)logT
    stats1 = @timed num = interpolate!(Ni, ωs[1:2*Nt], Nys[1:2*Nt])
    stats2 = @timed den = interpolate!(Di, ωs[1:2*Dt], Dys[1:2*Dt])
    push!(_runtime_vdhl_benortiwari, stats1.time)
    push!(_runtime_vdhl_benortiwari, stats2.time)
    # backward dilation,
    # substitute (x0,x1,x2,...xn) = (inv(d0)x0,inv(d1)x1,...,inv(dn)xn)
    # n*log(q), and T
    undilated = gens(Rhom) .* map(inv, dilation) 
    stats1 = @timed num = evaluate(num, undilated)
    stats2 = @timed den = evaluate(den, undilated)
    push!(_runtime_vdhl_pointsmanip, stats1.time)
    push!(_runtime_vdhl_pointsmanip, stats2.time)
    # dehomogenization,
    # substitute (x0,x1,x2...,xn) = (1,x1,x2...,xn),
    xs0 = [one(R), gens(R)...]
    num = evaluate(num, xs0)
    den = evaluate(den, xs0)
    # normalize by the trailing_coefficient,
    # T*n*log(q)
    normalization_factor = trailing_coefficient(den)
    num = map_coefficients(c -> div(c, normalization_factor), num)
    den = map_coefficients(c -> div(c, normalization_factor), den)
    num, den
end
