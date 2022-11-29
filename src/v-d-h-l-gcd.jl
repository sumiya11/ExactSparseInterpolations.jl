#=
    The implementation of the GCD algorithm from
    "On sparse interpolation of rational functions and gcds".
    https://dl.acm.org/doi/10.1145/3466895.3466896
=#

# vanDerHoevenLecerfGCD - an object used for GCD computation.
# ring - polynomial ring K[x1...xn],
# N,D - bounds on the degrees of numerator and denominator,
# univariate_poly_interpolator - an object that can interpolate dense polynomials,
# multivariate_poly_interpolator - an object that can interpolate polynomials;
mutable struct vanDerHoevenLecerfGCD{Ring, I1<:AbstractPolynomialInterpolator, I2<:AbstractPolynomialInterpolator}
    ring::Ring
    N::Int
    D::Int
    univariate_poly_interpolator::I1
    multivariate_poly_interpolator::I2
end

# Creates a vanDerHoevenLecerfGCD object from the given polynomial ring 
# and the given degree bounds
function vanDerHoevenLecerfGCD(ring, N::Integer, D::Integer; 
        univariate_poly_interpolator=Newton(univariatize(Nemo.PolyRing, ring), d=max(N, D)),
        multivariate_poly_interpolator=BenOrTiwari(homogenize(ring)))
    @assert N >= 0 && D >= 0
    vanDerHoevenLecerfGCD(ring, N, D, univariate_poly_interpolator, multivariate_poly_interpolator)
end

function gcd!(vdhlg::vanDerHoevenLecerfGCD, P, Q)
    R = vdhlg.ring
    K = base_ring(R)
    N, D = vdhlg.N, vdhlg.D
    upi = vdhlg.univariate_poly_interpolator
    P_interpolator = copy(upi)
    Q_interpolator = copy(upi)
    mpi = vdhlg.multivariate_poly_interpolator
    # polynomial ring k[x0,x1,..xn]
    Rhom = mpi.ring
    f = one(Rhom)
    # random variable shift
    shift = random_point(Rhom)
    i = 0
    success = false
    while !success
        # next point for sparse polynomial interpolation,
        # e.g., ωi = [2^i, 3^i, 5^i] in case BenOr-Tiwari algorithm is used
        ωi = next_point!(mpi)
        ω0 = first(ωi)
        i += 1
        # random points for dense polyomial interpolation,
        # f(ξ*x0,ξ*x1,ξ*x2,..., ξ*xn) for ξ in ξij
        ξij = distinct_points(K, max(N, D) + 1)
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
        Pij = map(P, ωξsijh)
        Qij = map(Q, ωξsijh)
        # multiply by the correcton factor
        Pij = map(cξ -> cξ[1] * (ω0*cξ[2] + shift[1])^(N), zip(Pij, ξij))
        Qij = map(cξ -> cξ[1] * (ω0*cξ[2] + shift[1])^(D), zip(Qij, ξij))
        # interpolate both polynomials densely in K[ξ]
        p = interpolate!(P_interpolator, ξij, Pij)
        q = interpolate!(Q_interpolator, ξij, Qij)
        empty!(P_interpolator)
        empty!(Q_interpolator)
        # compute the gcd of P and Q in K[ξ].
        g = Nemo.gcd(p, q)
        # normalize the gcd
        g = divexact(g, trailing_coefficient(g))
        lcg = leading_coefficient(g)
        success, f = next!(mpi, next_point!(mpi), lcg)
        if i > 2^10
            throw(ErrorException("Something bad happened in vanDerHoevenLecerfGCD"))
        end
    end
    xs0 = gens(R)
    xs0 = [one(R), xs0...]
    f = evaluate(f, xs0)
    f = divexact(f, leading_coefficient(f))
    f
end
