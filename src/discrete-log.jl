# Discrete logarithms in Z/Zp

#=
    Solves for xi:
    a1^x1*a2^x2*...*a3^x3 = y (mod p)
    
    ai are some (distinct) generators of the multiplicative subgroup
    of the finite field of prime characteristic p.

    O(n*log(p)*DL(p)),
    where DL(p) is the time complexity of the discrete log algorithm
    and n is the number of variables.
=#
function multibase_discrete_log(ai::Vector{I}, y::I, p::T) where {I, T}
    # solve for mi, variable by variable, the following:
    #   
    #   a1^m1*(a2a3...a4) = 1 (mod p)
    #   a2^m2*(a3...a4)   = 1 (mod p)
    #   ...
    #   an^mn             = y (mod p)
    # 
    # then xi = mi + i - 1
    n = length(ai)
    # vector of integers, the memory is not initilized
    xi = Vector{Int}(undef, n)
    # Z is inv(a2)*inv(a3)*...*inv(an)
    Z = prod(inv(ai[i]) for i in 2:n)
    @info "" Z
    for i in 1:n-1
        # ai^m = Z (mod p) for some m
        m = discrete_log(ai[i], Z, p)
        @info "" ai[1] Z p m
        @assert ai[i]^m == Z
        Z = Z*ai[i + 1]
        xi[i] = m + i - 1
    end
    xi[n] = n - 1 + discrete_log(ai[n], y, p)
    xi
end

# currently this is a random number 
_threshold_direct() = 2^4

# Solves a^x = y (mod p) for x
# ord is the order of a in Z/Zp.
function discrete_log(a::I, y::I, ord::T, factors) where {I, T}
    if ord < _threshold_direct_discreet_log()
        direct_discrete_log(a, y)
    else
        if ord == Nemo.ord(parent(a))
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
    m = ceil(Int, sqrt(ord))
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

function distinct_generators(k, n)
    ch = Int(characteristic(k))
    ord = ch - 1
    factors = Primes.factor(Vector, ord)
    v = zeros(k, n)
    jj = 2
    for i in 1:n
        while jj < ch
            # kj = k(jj)
            kj = rand(k)
            if generates_mult_group(ord, factors, kj)
                v[i] = kj
                jj += 1
                break
            end
            jj += 1
        end
        jj >= ch && error("The characteristic of the base field is too small, sorry")
    end
    v
end

# Returns alphanp2, a generator of K, such that
# alphanp2/alphanp1 is itself a generator of K,
# and alphanp2 is not in alphan
function distinct_generator(K, alphan, alphanp1)
    ord = Int(order(K) - 1)
    factors = Primes.factor(Vector, ord)
    i = 0
    while true
        i += 1
        i > ord && error("Could not find a distinct generator")
        gamma = rand(K)
        alphanp2 = gamma*alphanp1 
        alphanp2 in alphan && continue
        if generates_mult_group(ord, factors, gamma)
            return alphanp2
        end
    end
    return one(K)
end

function generates_mult_group(ord, factors, kj)
    iszero(kj) && return false
    for p in factors
        d = div(ord, p)
        if isone(kj^d)
            return false
        end
    end
    true
end

# Returns a primitive r-th root of unity in K.
function root_of_unity(K, r)
    ord = order(K)
    factors = Primes.factor(Vector, ord)
    while true
        a = rand(K)
        if generates_mult_group(ord, factors, a)
            return a^div(ord, r)
        end
    end
    # bad branch
    return one(K)
end

# Returns a primitive root of unity in K
# of order >= r and the order itself
function approx_root_of_unity(K, r)
    ord = Int(order(K)) - 1
    @assert r <= ord
    g = first(distinct_generators(K, 1))
    r == ord && return (ord, g)
    factors = Primes.factor(Vector, ord)
    i = 1
    ord1 = ord
    while ord1 >= r && i <= length(factors)
        @assert ord1 >= r
        ord1 = div(ord1, factors[i])
        i += 1
    end
    ord1 = ord1*factors[i - 1]
    @assert ord1 >= r
    (ord1, g^div(ord, ord1))
end
