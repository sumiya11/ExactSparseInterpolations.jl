
mutable struct PolyQQInterpolator{Ring}
    # univariate polynomial ring
    ring::Ring
    # the number of terms in the interpolant
    T::Int
    # the degree of the interpolant
    D::Int
    function PolyQQInterpolator(ring::Ring, T::Integer, D::Integer) where {Ring}
        @assert T >= 0 && D >= 0
        @assert nvars(ring) == 1
        @assert base_ring(ring) == Nemo.QQ "Base ring must be Nemo.QQ"
        new{Ring}(ring, T, D)
    end
end

function find_reduction_moduli(T, D)
    moduli = Vector{Int}()
    # p = 2T + 1
    p = 2^20
    product = fmpz(1)
    while product < fmpz(D)
        p = Primes.nextprime(p + 1)
        product *= p
        push!(moduli, p)
    end
    moduli
end

function interpolate!(I::PolyQQInterpolator, blackbox)
    moduli = find_reduction_moduli(I.T, I.D)
    imodp = []
    # @info "" moduli
    for p in moduli
        FF = Nemo.GF(p)
        Rp, _ = PolynomialRing(FF, ["xp"])
        bot = FasterBenOrTiwari(Rp, I.T)
        ω = randomgenerator(FF)
        ωi = [[ω^i] for i in 0:2I.T-1]
        bbp = change_base_ring(FF, blackbox)
        fi = map(bbp, ωi)
        fp = interpolate!(bot, ωi, fi)
        push!(imodp, fp)
    end
    # @info "" imodp
    x = gens(I.ring)[1]
    fstar = zero(I.ring)
    Ct = map(f -> strictly_unique(term -> lift(coeff(term, 1)), collect(terms(f))), imodp)
    Ck = map(f -> map(term -> lift(coeff(term, 1)), f), Ct)
    Call = union(Ck...)
    for (i, c) in enumerate(Call)
        inds = filter(j -> c in Ck[j], 1:length(Ck))
        totaldeg = lcm(map(j -> moduli[j] - 1, inds))
        # @info "" i c
        # @info "" inds totaldeg
        if totaldeg > I.D
            emodr = map(j -> degree(Ct[j][findfirst(term -> lift(coeff(term, 1)) == c, Ct[j])], 1), inds)
            modr = map(j -> moduli[j] - 1, inds)
            # @info "" emodr modr
            e = Nemo.crt(emodr, modr)
            newmonom = c*x^e
            # @warn "Success!" newmonom 
            fstar += newmonom
        end
    end
    fstar
end
