
# Computes the discrete Fourier transform of vector xs, 
# using the root of unity ω of order r, that is,
#   yi = sum_j xs[j] ω^(i*j)
function dft(xs::Vector{T}, ω::T, r::Integer) where {T}
    n = length(xs)
    @assert n == r
    K = parent(ω)
    ys = zeros(K, n)
    for i in 0:n-1
        for j in 0:n-1
            ys[i + 1] += xs[j + 1] * ω^(i*j)
        end
    end
    ys
end

# Computes the inverse discrete Fourier transform of vector xs,
# using the root of unity ω of order r, that is,
#   xi = -C * sum_j ys[j] ω^(-i*j)
function invdft(ys::Vector{T}, ω::T, r::Integer) where {T}
    n = length(ys)
    K = parent(ω)
    @assert isone(ω^r)
    C = div(order(K) - 1, r)
    xs = zeros(K, n)
    for i in 0:n-1
        for j in 0:n-1
            xs[j + 1] += -C * ys[i + 1] * ω^(-i*j)
        end
    end
    xs
end

function fftinterpolate(R, blackbox, r)
    K = base_ring(R)
    ord, ω = approx_root_of_unity(K, r)
    @info "" r ord ω
    ωi = [ω^i for i in 0:ord-1]
    @assert isone(ωi[1]) && !any(isone, ωi[2:end])
    fi = map(blackbox, ωi)
    ci = invdft(fi, ω, ord)
    R(ci)
end

# 1048583

function select_moduli(K, T, D)
    rs = Vector{Any}()
    product = fmpz(1)
    if K == Nemo.GF(73)
        factors = [8, 9]
        generator = K(26)
        candidates = [
            (generator^(8), 9),
            (generator^(9), 8),
        ]
        i = 1
        while i <= length(candidates) && product < D
            ω, ord = candidates[i]
            product *= fmpz(ord)
            @info "" ord ω product
            push!(rs, (ω, ord))
            i += 1
        end
        product >= D && return rs
        error("Too small field")
    end
    if K == Nemo.GF(239)
        factors = [2, 7, 17]
        generator = K(37)
        candidates = [
            (generator^(2*7), 17),
            (generator^(17), 14),
        ]
        i = 1
        while i <= length(candidates) && product < D
            ω, ord = candidates[i]
            product *= fmpz(ord)
            @info "" ord ω product
            push!(rs, (ω, ord))
            i += 1
        end
        product >= D && return rs
        error("Too small field")
    end
    while product < D
        ord, ω = approx_root_of_unity(K, T)
        product *= fmpz(ord)
        @info "" ord ω product
        push!(rs, (ω, ord))
    end
    for i in 1:length(rs)
        for j in i+1:length(rs)
            if !isone(gcd(rs[i][2], rs[j][2]))
                return select_moduli(K, T, D)
            end
        end
    end
    rs
end

function evaluate_in_cyclic_ext_1(R, blackbox, ω, ord)
    x = gens(R)
    modulo = x[1]^ord - 1
    mod(blackbox(x), modulo)
end

function evaluate_in_cyclic_ext_2(R, blackbox, ω, ord)
    ωi = [[ω^i] for i in 0:ord-1]
    @assert isone(ωi[1][1]) # && !any(isone, ωi[2:end])
    fi = map(blackbox, ωi)
    ci = invdft(fi, ω, ord)
    ms = map(i -> [i], filter(i -> !iszero(ci[i + 1]), 0:length(ci)-1))
    ci = map(i -> ci[i[1] + 1], ms)
    R(ci, ms)
end

function strictly_unique(f, arr)
    seen = Dict()
    for x in arr
        if haskey(seen, f(x))
            delete!(seen, f(x))
        else
            seen[f(x)] = x
        end
    end
    collect(values(seen))
end

function approximate_interpolate(R, T, D, blackbox)
    K = base_ring(R)
    x = gens(R)[1]
    @assert order(K) > 2T
    # select rk, such that rk > T, prod(rk) > D,
    # order(rk) >~ T
    rs = select_moduli(K, T, D)
    @info "" rs
    frs = map(
        rk -> evaluate_in_cyclic_ext_2(R, blackbox, rk[1], rk[2]),
        rs
    )
    @info "" frs
    fstar = zero(R)
    Ct = map(f -> strictly_unique(term -> coeff(term, 1), collect(terms(f))), frs)
    Ck = map(f -> map(term -> coeff(term, 1), f), Ct)
    Call = union(Ck...)
    @info "" Ct Ck Call
    for (i, c) in enumerate(Call)
        inds = filter(j -> c in Ck[j], 1:length(Ck))
        totaldeg = lcm(map(j -> rs[j][2], inds))
        @info "" i c
        @info "" inds totaldeg
        if totaldeg > D
            emodr = map(j -> degree(Ct[j][findfirst(term -> coeff(term, 1) == c, Ct[j])], 1), inds)
            modr = map(j -> rs[j][2], inds)
            @info "" emodr modr
            e = Nemo.crt(emodr, modr)
            newmonom = c*x^e
            @warn "Success!" newmonom 
            fstar += newmonom
        end
    end
    fstar
end
