#=
    The implementation of the algorithm from
    "On sparse interpolation of rational functions and gcds".
    https://dl.acm.org/doi/10.1145/3466895.3466896
=#

# vanDerHoevenLecerf - an object used for interpolation.
# ring - polynomial ring K[x1...xn],
# N,D - bounds on the degrees of numerator and denominator,
# univariate_rational_interpolator - an object that can interpolate rational functions,
# multivariate_poly_interpolator - an object that can interpolate polynomials;
mutable struct vanDerHoevenLecerf{Ring, I1<:AbstractRationalInterpolator, I2<:AbstractPolynomialInterpolator} <: AbstractRationalInterpolator
    ring::Ring
    N::Int
    D::Int
    univariate_rational_interpolator::I1
    multivariate_poly_interpolator::I2
end

# Creates a vanDerHoevenLecerf object from the given polynomial ring 
# and the given degree bounds
function vanDerHoevenLecerf(ring, N::Integer, D::Integer; 
        univariate_rational_interpolator=Cauchy(univariatize(Nemo.PolyRing, ring), N, D),
        multivariate_poly_interpolator=BenOrTiwari(homogenize(ring)))
    @assert N >= 0 && D >= 0
    vanDerHoevenLecerf(ring, N, D, univariate_rational_interpolator, multivariate_poly_interpolator)
end

# Given a vanDerHoevenLecerf object and a blackbox function,
# returns a tuple of polynomials (P, Q), 
# such that P/Q interpolates the blackbox function
function interpolate!(vdhl::vanDerHoevenLecerf, blackbox)
    R = vdhl.ring
    xs = Nemo.gens(R)
    K = base_ring(R)
    N, D = vdhl.N, vdhl.D
    uri = vdhl.univariate_rational_interpolator
    mpi = vdhl.multivariate_poly_interpolator
    # evaluated leading coefficients of numerator and denominator
    P_coeff = Vector{elem_type(K)}(undef, 0)
    Q_coeff = Vector{elem_type(K)}(undef, 0)
    # interpolated leading coefficients of numerator and denominator
    P_interpolated = Ref{elem_type(R)}(zero(R))
    Q_interpolated = Ref{elem_type(R)}(zero(R))
    # corresponding interpolators
    P_interpolator = copy(mpi)
    Q_interpolator = copy(mpi)
    # polynomial ring k[x0,x1,..xn]
    Rhom = mpi.ring
    # random variable shift
    shift = random_point(Rhom)
    # multivariate polynomial interpolation points
    ωs = Vector{Vector{elem_type(K)}}(undef, 0)
    i = 0
    all_interpolated = false
    while !all_interpolated
        # next point for sparse polynomial interpolation,
        # e.g., ωi = [2^i, 3^i, 5^i] in case BenOr-Tiwari algorithm is used
        ωi = next_point!(mpi, increment=true)
        push!(ωs, ωi)
        i += 1
        # homogenizing component ω0
        ω0 = ωi[1]
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
        # evaluate the blackbox
        fij = map(blackbox, ωξsijh)
        # multiply by the correcton factor ω0^(N - D)
        fij = map(cξ -> cξ[1] * (ω0*cξ[2] + shift[1])^(N - D), zip(fij, ξij))
        # interpolate the numerator and the denominator densely.
        P, Q = interpolate!(uri, ξij, fij)
        @assert isone(trailing_coefficient(Q))
        empty!(uri)
        # n and d are the true dense degrees of the numerator and denominator
        n, d = degree(P), degree(Q)
        # store coefficients of dense interpolation of P and Q
        # in P_coeff and Q_coeff respectively
        for (poly, cfs) in ((P, P_coeff), (Q, Q_coeff))
            push!(cfs, leading_coefficient(poly))
        end
        # interpolate the leading coefficient 
        # simultaneously in numerator and denominator
        all_interpolated = true
        for (cfs, interpolated, interpolator) in (
            (P_coeff, P_interpolated, P_interpolator),
            (Q_coeff, Q_interpolated, Q_interpolator)
        )   
            y_point = last(cfs)
            success, f = next!(interpolator, next_point!(interpolator), y_point)
            interpolated[] = f
            all_interpolated = all_interpolated && success
        end
        if i > 2^20
            throw(ErrorException("Something bad happened in vanDerHoevenLecerf"))
        end
    end
    # P and Q are the interpolated leading coefficients in K[x0,x1,...xn]
    P = P_interpolated[]
    Q = Q_interpolated[]
    # dehomogenize P and Q
    xs0 = gens(R)
    xs0 = [one(R), xs0...]
    P = evaluate(P, xs0)
    Q = evaluate(Q, xs0)
    normalization_factor = trailing_coefficient(Q)
    P = map_coefficients(c -> div(c, normalization_factor), P)
    Q = map_coefficients(c -> div(c, normalization_factor), Q)
    P, Q
end
