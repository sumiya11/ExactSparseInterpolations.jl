#=
    Berlekamp-Massey algorithm.

    Notation in the implementation is taken from Algorithm 1 from
        https://arxiv.org/pdf/2004.01463.pdf
=#

mutable struct BerlekampMassey{Ring, Point, Poly} <: AbstractPolynomialInterpolator
    ring::Ring
    a::Vector{Point}
    Δ::Point
    Λ::Poly
    B::Poly
    L::Int
    function BerlekampMassey{Ring, Point, Poly}(ring::Ring) where {Ring, Point, Poly}
        new{Ring, Point, Poly}(
            ring, zeros(base_ring(ring), 0),
            one(base_ring(ring)), one(ring), one(ring), zero(Int)
        ) 
    end
    function BerlekampMassey{Ring, Point, Poly}(ring::Ring, a::Vector{Point}, Δ::Point, Λ::Poly, B::Poly, L::Int) where {Ring, Point, Poly}
        new{Ring, Point, Poly}(ring, a, Δ, Λ, B, L) 
    end
end

BerlekampMassey(ring::Ring) where {Ring} = BerlekampMassey{Ring, elem_type(base_ring(ring)), elem_type(ring)}(ring)

function Base.copy(bm::BerlekampMassey{Ring, Point, Poly}) where {Ring, Point, Poly}
    BerlekampMassey{Ring, Point, Poly}(bm.ring, copy(bm.a), bm.Δ, bm.Λ, bm.B, bm.L)
end

function Base.empty!(bm::BerlekampMassey)
    empty!(bm.a)
    bm.Λ = one(bm.ring)
    bm.B = one(bm.ring)
    bm.L = zero(Int)
    bm.Δ = one(base_ring(bm.ring))
    bm
end

# Given a BerlekampMassey algorithm instance and an interpolation point
# returns a tuple (flag, poly),
# where flag is a boolean variable set to true if the discrepancy is zero,
# and poly is the generating polynomial.
function next!(bm::BerlekampMassey{Ring, Point, Poly}, x::Point) where {Ring, Point, Poly}
    # z is the generator of the univariate polynomial ring
    z = gen(bm.ring)
    # K is the ground field
    K = base_ring(bm.ring)
    push!(bm.a, x)
    r = length(bm.a)
    # let bm.Λ be the generating polynomial of the recurrent sequence computed so far;
    # then λs is [bm.Λ_1, bm.Λ_2, ..., bm.Λ_d, 0, ..., 0],
    # where the number of zeroes is r - deg(λs) - 1
    λs = append!(collect(coefficients(bm.Λ)), zeros(K, r-degree(bm.Λ)-1))
    # discrepancy Δr
    Δr = sum(map(prod, zip(reverse(λs), bm.a)))
    if iszero(Δr)
        bm.B = z*bm.B
    elseif 2*bm.L < r
        bmB = bm.B
        bm.B = bm.Λ
        bm.Λ = bm.Λ - Δr//bm.Δ*z*bmB
        bm.L = r - bm.L
        bm.Δ = Δr
    else
        bm.Λ = bm.Λ - Δr//bm.Δ*z*bm.B
        bm.B = z*bm.B
    end
    iszero(Δr), bm.Λ
end
