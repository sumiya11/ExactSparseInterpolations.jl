
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
    ωi = [ω^i for i in 0:ord-1]
    @assert isone(ωi[1]) && !any(isone, ωi[2:end])
    fi = map(blackbox, ωi)
    ci = invdft(fi, ω, ord)
    R(ci)
end
