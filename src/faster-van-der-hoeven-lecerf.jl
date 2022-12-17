
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
    # interpolator for univariate functions
    cauchy::FasterCauchy{UnivRing}
    # polynomial interpolators for the LCs of numerator/denominator
    Ni::FasterMultivariateBenOrTiwari{Ring}
    Di::FasterMultivariateBenOrTiwari{Ring}

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
        Ni = FasterMultivariateBenOrTiwari(ringhom, Nt, Nds)
        Di = FasterMultivariateBenOrTiwari(ringhom, Dt, Dds)
        new{Ring,typeof(Runiv)}(
            ring, Nd, Dd, Nds, Dds, Nt, Dt, cauchy, Ni, Di
        )
    end
end

function interpolate!(vdhl::FasterVanDerHoevenLecerf, blackbox)
    R = vdhl.ring
    K = base_ring(R)
    Nd, Dd = vdhl.Nd, vdhl.Dd
    Nt, Dt = vdhl.Nt, vdhl.Dt
    Ni, Di = vdhl.Ni, vdhl.Di
    cauchy = vdhl.cauchy
    T = max(Nt, Dt)
    #
    Rhom = Ni.ring
    shift = random_point(Rhom)
    #
    ω = startingpoint(Ni)
    # @info "" ω shift Nt Dt T
    # 
    ωs = map(i -> ω .^ i, 0:2T-1)
    Nys = Vector{elem_type(K)}(undef, 2T)
    Dys = Vector{elem_type(K)}(undef, 2T)
    # @warn "" Ni.Ds Di.Ds Ni.Dsubs Di.Dsubs
    for i in 0:2T-1
        ωi = ωs[i + 1]
        # @info "" i Nd + Dd + 2 ωi
        ξij = distinct_points(K, Nd + Dd + 2)
        # "substitute" ωi,
        # f(ξ*ωi0, ξ*ωi1,ξ*ωi2,..., ξ*ωin) for ξ in ξij
        ωξij = [ωi .* ξ for ξ in ξij]
        # shift in each variable,
        # f(ξ*ωi0 + s0, ξ*ωi1 + s1,ξ*ωi2 + s2,..., ξ*ωin + sn) for ξ in ξij
        ωξsij = [ωξ .+ shift for ωξ in ωξij]
        # homogenize
        # ωξsijh = { (ξ*ωi1 + s1)/(ξ*ωi0 + s0), (ξ*ωi2 + s2)/(ξ*ωi0 + s0)... } for ξ in ξij
        ωξsijh = [(ωξs .// ωξs[1])[2:end] for ωξs in ωξsij]
        # evaluate the blackbox
        fij = map(blackbox, ωξsijh)
        # multiply by the correction factor ω0^(N - D)
        # fij = map(cξ -> cξ[1] * (ω0*cξ[2] + shift[1])^(Nd - Dd), zip(fij, ξij))
        fij = map(cξ -> cξ[1] * cξ[2][1]^(Nd - Dd), zip(fij, ωξsij))
        # interpolate the numerator and the denominator densely.
        N, D = interpolate!(cauchy, ξij, fij)
        # @info "" N D
        @assert isone(trailing_coefficient(D))
        Nys[i + 1] = leading_coefficient(N)
        Dys[i + 1] = leading_coefficient(D)
    end
    #
    # @info "" ωs[1:2*Dt] Dys[1:2*Dt]
    num = interpolate!(Ni, ωs[1:2*Nt], Nys[1:2*Nt])
    den = interpolate!(Di, ωs[1:2*Dt], Dys[1:2*Dt])
    # @info "" num den
    #
    xs0 = [one(R), gens(R)...]
    num = evaluate(num, xs0)
    den = evaluate(den, xs0)
    # 
    normalization_factor = trailing_coefficient(den)
    num = map_coefficients(c -> div(c, normalization_factor), num)
    den = map_coefficients(c -> div(c, normalization_factor), den)
    num, den
end
