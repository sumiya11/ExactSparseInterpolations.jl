#=
    Interpolation variable by variable [Zippel, 1990].

    Multivariate, Dense-Sparse, No-Early-Termination, Blackbox.

    The algorithm comes in two flavors:
    MC Probabilistic, and Deterministic

    Probabilistic version:
        O(d*n*t) queries to blackbox (exactly, d*n*t - v*t - d*t ?)
        O(d*n*t^2) runtime 
    Deterministic version:

=#

mutable struct Zippel{Ring, Point, T3} <: AbstractPolynomialInterpolator
    ring::T1
    d::Int
    t::Int
    T::Int
    start_point::Vector{T2}
    dense_univariate_interpolator::T3
end

function Zippel(ring, d::Integer;
                t::Integer=-1, T::Integer=-1,
                start_point=random_point(ring),
                dense_univariate_interpolator=Newton(univariatize(ring), d=d))
    Zippel(ring, d, t, T, start_point, dense_univariate_interpolator)
end

function skeleton_poly(poly)
    collect(exponent_vectors(poly))
end

function skeleton_eval(poly, p)
    [prod((p .^ x)) for x in poly]
end

function interpolate!(z::Zippel{T1, T2, T3}, blackbox) where {T1, T2, T3}
    R = parent(z)
    v = nvars(R)
    K = base_ring(R)
    sp = z.start_point
    d = z.d
    t, T = z.t, z.T
    dui = z.dense_univariate_interpolator
    sp1 = sp[1]
    xs = Vector{elem_type(K)}(undef, d)
    fx = Vector{elem_type(K)}(undef, d)
    for i in 1:d
        xs[i] = generic_point(K)
        sp[1] = xs[i]
        fx[i] = blackbox(sp)
    end
    sp[1] = sp1
    P1 = interpolate(dui, xs, fx)
    S = skeleton_poly(P1)
    # @info "" P1 S
    # @info "" v d
    for i in 2:v
        C = []
        push!(C, collect(coefficients(P1)))
        R = zeros(elem_type(K), d)
        R[1] = sp[i]
        for j in 2:d
            rj = generic_point(K)
            R[j] = rj
            t = length(S)
            Lspace = MatrixSpace(K, t, t)
            L = zero(Lspace)
            b = zeros(elem_type(K), t)
            for k in 1:t
                lambda_k = [generic_point(K) for _ in 1:i-1]
                Slambda = skeleton_eval(S, lambda_k)
                # @info "generated point" lambda_k Slambda
                for kk in 1:t
                    L[k, kk] = Slambda[kk]
                end
                b[k] = blackbox([lambda_k..., rj, sp[i+1:end]...])
            end
            # @info "" j L b
            Linv = Nemo.inv(L)
            push!(C, Linv*b)
        end
        # @info "" i R C
        Ys = []
        for ii in 1:length(C[1])
            Cs = [c[ii] for c in C]
            Yc = interpolate(dui, R, Cs)
            push!(Ys, Yc)
            # @info "" R Cs Yc
        end
        Ri, _ = PolynomialRing(K, map(string, gens(parent(z)))[1:i])
        P1 = Ri(collect(coefficients(P1)), map(x -> [x..., 0], collect(exponent_vectors(P1))))
        Ys = map(y -> Ri(collect(coefficients(y)), map(x -> [zeros(Int, i-1)..., x...], collect(exponent_vectors(y)))), Ys)
        P1 = sum(collect(monomials(P1)) .* Ys)
        S = skeleton_poly(P1)
        # @warn "" P1 S
    end
    change_base_ring(K, P1, parent=parent(z))
end
