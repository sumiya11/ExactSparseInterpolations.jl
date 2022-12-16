# Discrete logarithms in F_q^n

# Stores the field K and some information about it
# to speed up discrete logarithms
mutable struct PrecomputedField{Field, T, I}
    K::Field
    ordmult::T
    factors::Vector{Pair{I, Int}}
    extensiondeg::Int
end

function PrecomputedField(K::Field) where {Field}
    ordmult = Int(Nemo.order(K) - 1)
    factors = collect(Primes.factor(Dict, ordmult))
    PrecomputedField(K, ordmult, factors, Nemo.degree(K))
end

# currently this is a random number 
_threshold_direct() = 2^4

function discrete_log(a, y, K::PrecomputedField)
    [discrete_log(a[1], y, K.ordmult, K.factors)]
end

# Solves a^x = y (mod p) for x
# ord is the order of a in Z/Zp, factors is an array of factors of p-1
function discrete_log(a::I, y::I, ord::T, factors) where {I, T}
    if ord < _threshold_direct()
        direct_discrete_log(a, y)
    else
        if ord == Nemo.order(parent(a))-1
            pohlig_hellman_discrete_log(a, y, ord, factors)
        else
            babystep_giantstep_discrete_log(a, y, ord)
        end
    end
end

# Solves a^x = y (mod p) for x using the Pohlig-Hellman algorithm.
# 
# ord is the order of a in Z/Zp.
# factors is a dictionary of prime factors of ord (with multiplicities).
# xibuf and pibuf are buffers used to store intermediate results.
function pohlig_hellman_discrete_log(a::I, y::I, ord::T, factors;
        xibuf=Vector{T}(undef, length(factors)),
        pibuf=Vector{T}(undef, length(factors))) where {I, T}
    @inbounds for i in 1:length(factors)
        (pi, di) = factors[i]
        ai, yi = a, y
        xi = zero(T)
        cc = one(ai)
        for j in 0:di-1
            pij = pi^j
            cij = div(ord, pij*pi)
            aij = ai^(cij*pij)
            yij = (yi*inv(cc))^cij
            xij = discrete_log(aij, yij, pi, factors)
            tij = xij*pij
            cc *= ai^tij
            xi += tij
        end
        xibuf[i] = xi
        pibuf[i] = pi^di
    end
    Nemo.crt(xibuf, pibuf)
end

# Solves a^x = y (mod p) for x using the baby-step giant-step algorithm.
# ord is the order of a in Z/Zp.
function babystep_giantstep_discrete_log(a::I, y::I, ord::T) where {I, T}
    # the size of the giant step
    m = ceil(BigInt, sqrt(ord))
    # baby steps
    baby = Dict{I, T}()
    ai = one(a)
    for i in 0:m-1
        baby[ai] = i
        ai *= a
    end
    # find a match
    ainvm = inv(a)^m
    giantstep = y
    i = 0
    while !haskey(baby, giantstep)
        giantstep *= ainvm
        i += 1
    end
    i*m + baby[giantstep]
end

# Solves a^x = y (mod p).
function direct_discrete_log(a::I, y::I) where {I}
    i = 0
    ai = one(a)
    while ai != y
        ai *= a
        i += 1
    end
    i
end
