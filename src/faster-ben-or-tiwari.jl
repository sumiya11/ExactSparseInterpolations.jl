
# Measure running time of each step!
const _runtime_benortiwari_padé = Float64[]
const _runtime_benortiwari_roots = Float64[]
const _runtime_benortiwari_dlog = Float64[]
const _runtime_benortiwari_vandermonde = Float64[]

function _runtime_benortiwari_dump()
    ans = (
        padé=minimum(_runtime_benortiwari_padé),
        roots=minimum(_runtime_benortiwari_roots),
        dlog=minimum(_runtime_benortiwari_dlog),
        vandermonde=minimum(_runtime_benortiwari_vandermonde),
    )
    empty!(_runtime_benortiwari_padé)
    empty!(_runtime_benortiwari_roots)
    empty!(_runtime_benortiwari_dlog)
    empty!(_runtime_benortiwari_vandermonde)
    ans
end

# Interpolate polynomials in O(M(n)logn),
# where M(n) is the complexity of polynomial multiplication
mutable struct FasterBenOrTiwari{Ring} <: AbstractPolynomialInterpolator
    # univariate polynomial ring
    ring::Ring
    # the number of terms in the interpolant
    T::Int

    function FasterBenOrTiwari(ring::Ring, T::Integer) where {Ring}
        @assert T >= 0
        @assert nvars(ring) == 1
        new{Ring}(ring, T)
    end
end

# Same as above, but multivariate.
# Currently uses Kronecker substitutions.
mutable struct FasterMultivariateBenOrTiwari{Ring} <: AbstractPolynomialInterpolator
    # multivariate polynomial ring
    ring::Ring
    # the number of terms in the interpolant
    T::Int
    # the vector of partial degrees of the interpolant
    Ds::Vector{Int}
    # the vector of degrees used in Kronecker substitution
    Dsubs::Vector{BigInt}
    function FasterMultivariateBenOrTiwari(ring::Ring, T::Integer, Ds::Vector{<:Integer}) where {Ring}
        @assert T >= 0
        @assert all(>=(0), Ds)
        @assert length(Ds) == nvars(ring)
        Dsubs = subsdegrees(Ds)
        K = base_ring(ring)
        Dsubs[end] >= order(K) && @warn "In Kronecker substitution the field order is too small" Dsubs order(K)
        new{Ring}(ring, T, Ds, Dsubs)
    end
end

const AbstractBenOrTiwari = Union{FasterBenOrTiwari,FasterMultivariateBenOrTiwari}

# Returns [1, D1, D1D2, ..., D1...Dn] to be used
# in Kronecker substitution.
# (probably, adding 1 to each entry)
function subsdegrees(Ds::Vector{<:Integer})
    ans = Vector{BigInt}(undef, length(Ds) + 1)
    ans[1] = one(BigInt)
    for i in 1:length(Ds)
        ans[i+1] = ans[i]*(Ds[i] + 1)
    end
    ans
end

subsbackward(fbot::FasterBenOrTiwari, monoms) = map(e -> map(Int, e), monoms)
subsbackward(fmbot::FasterMultivariateBenOrTiwari, monoms::Vector{Vector{I}}) where {I} = map(m -> subsbackward(fmbot, m), monoms)
function subsbackward(fmbot::FasterMultivariateBenOrTiwari, monom::Vector{I}) where {I}
    @assert length(monom) == 1
    m = monom[1]
    Dsubs = fmbot.Dsubs
    n = length(Dsubs) - 1
    ans = Vector{Int}(undef, n)
    for i in n:-1:1
        ans[i] = div(m, Dsubs[i])
        m = m - ans[i]*Dsubs[i]
    end
    ans
end

# Returns an appropriate starting point 
# for the geometric sequence ω^0, ω^1, ω^2, ...
function startingpoint(fbot::FasterBenOrTiwari)
    R = fbot.ring
    K = base_ring(R)
    n = nvars(R)
    @assert n == 1
    # if univariate
    if degree(K) == 1
        # if the field is actually Z/Zp,
        # return an element of order p-1
        [randomgenerator(K)]
    else
        # otherwise, return a primitive element
        [gen(K)]
    end
end

function startingpoint(fmbot::FasterMultivariateBenOrTiwari)
    # Do Kronecker substitution
    Dsubs = fmbot.Dsubs
    K = base_ring(fmbot.ring)
    g = randomgenerator(K)
    map(i -> g^Dsubs[i], 1:length(Dsubs)-1)
end

function interpolate!(fbot::AbstractBenOrTiwari, blackbox)
    T = fbot.T
    # the base of the geometric sequence ω^1, ω^2, ...
    ω = startingpoint(fbot)
    # generate the sequence
    # O(TlogT)
    ωs = map(i -> ω .^ i, 0:2T-1)
    # evaluate the blackbox function at the sequence
    # O(LT), where L is the cost of evaluating the blackbox function
    ys = map(blackbox, ωs)
    interpolate!(fbot, ωs, ys)
end

function interpolate!(fbot::AbstractBenOrTiwari, xs, ys)
    # check that the first degree is 0
    @assert all(isone, xs[1])
    @assert length(xs) == length(ys) == 2*fbot.T
    ω = xs[2]
    Rx = fbot.ring
    K = base_ring(Rx)
    T = fbot.T
    # construct the polynomial ys[1]z^1 + ys[2]z^2 + ... + ys[2T]z^(2T)
    # O(1)
    Rz, z = K["z"]
    sequence = Rz(ys)
    # find A/B such that A/B = sequence mod z^(2T) and degree(A) < T
    # O(M(T)logT)
    stats = @timed _, B, _ = Padé(sequence, z^(2T), T - 1)
    push!(_runtime_benortiwari_padé, stats.time)
    # @info "" ω B
    @assert degree(B) == T
    # assuming this is O(T logT^k logq^m) for some k and m, 
    # where q is the order of the base field
    stats = @timed mi = map(inv, Nemo.roots(B))
    push!(_runtime_benortiwari_roots, stats.time)
    @assert length(mi) == T
    # @info "" mi
    # find the monomials of the interpolant,
    # O(TlogTlogq), where q is the order of the base field
    # (note that this cost covers the case where K is not a prime field)
    # (assuming ord is smooth)
    PF = PrecomputedField(K)
    stats = @timed monoms = map(m -> discrete_log(ω, m, PF), mi)
    push!(_runtime_benortiwari_dlog, stats.time)
    # @info "" monoms
    # find the coefficients of the interpolant
    # by solving a T×T Vandermonde system
    # O(M(T)logT)
    stats = @timed coeffs = solve_transposed_vandermonde(Rz, view(mi, 1:T), view(ys, 1:T))
    push!(_runtime_benortiwari_vandermonde, stats.time)
    Rx(coeffs, subsbackward(fbot, monoms))
end
