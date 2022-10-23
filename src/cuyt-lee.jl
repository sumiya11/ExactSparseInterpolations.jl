
mutable struct CuytLee{T1, T2, T3} <: AbstractInterpolator
    ring::T1
    N::Int
    D::Int
    T::Int
    dense_univariate_interpolator::T2
    sparse_multivariate_interpolator::T3
end

function CuytLee(ring, N, D; 
        T=0, 
        dense_univariate_interpolator=RationalDifferenceNormalized(true_univariatize(ring), N, D),
        sparse_multivariate_interpolator=BenOrTiwari)
    CuytLee(ring, N, D, T, dense_univariate_interpolator, sparse_multivariate_interpolator)
end

function Nemo.interpolate(cl::CuytLee{T1}, blackbox) where {T1}
    R = cl.ring
    K = base_ring(R)
    N, D = cl.N, cl.D
    dui = cl.dense_univariate_interpolator
    smi = cl.sparse_multivariate_interpolator
    
    ω = geometric_point(R)
    all_interpolated = false
    Pints = [smi(R) for _ in 0:N]
    Qints = [smi(R) for _ in 1:D]
    Psuccess = [false for _ in 0:N]
    Qsuccess = [false for _ in 1:D]
    Pcoeffs = [one(R) for _ in 0:N]
    Qcoeffs = [one(R) for _ in 0:D]

    i = 1
    while !all_interpolated
        all_interpolated = true
        ξij = [generic_point(K) for _ in 0:N + D]
        ωi = ω .^ i
        ωξij = [ωi .* ξ for ξ in ξij]
        i += 1
        fij = map(blackbox, ωξij)
        P, Q = interpolate(dui, ξij, fij)
        @info "" P Q
        for (jj, (coeff, interpol)) in enumerate(zip(coefficients(P), Pints))
            success, f = next!(interpol, ωi, coeff)
            Pcoeffs[jj] = f
            Psuccess[jj] = success
            if !success
                all_interpolated = false
            end
        end
        @info Psuccess Pcoeffs
        for (jj, (coeff, interpol)) in enumerate(zip(coefficients(Q), Qints))
            success, f = next!(interpol, ωi, coeff)
            Qcoeffs[jj] = f
            Qsuccess[jj] = success
            if !success
                all_interpolated = false
            end
        end
        @info Qsuccess Qcoeffs
    end

    sum(Pcoeffs)-2, sum(Qcoeffs)-1
end
