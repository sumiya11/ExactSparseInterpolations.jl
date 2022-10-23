
mutable struct BerlekampMassey{I, R, Y}
    Ru::R
    a::Vector{I}
    Λ::Y
    B::Y
    L::Int
    Δ::I
    function BerlekampMassey{I, R, Y}(Ru::R; 
            a=zeros(I, 0)) where {I, R, Y}
        new{I, R, Y}(
            Ru, a,
            one(Ru), zero(Ru), zero(Int), one(base_ring(Ru))
        ) 
    end
end

function BerlekampMassey(Ru::R) where {R}
    br = base_ring(Ru)
    et = elem_type(Ru)
    BerlekampMassey{elem_type(br), R, et}(Ru)
end

function Base.empty!(bm::BerlekampMassey)
    empty!(bm.a)
    bm.Λ = one(bm.Ru)
    bm.B = zero(bm.Ru)
    bm.L = zero(Int)
    bm.Δ = one(base_ring(bm.Ru))
    bm
end

function next!(bm::BerlekampMassey{I, R, Y}, x::I) where {I, R, Y}
    z = gen(bm.Ru)
    push!(bm.a, x)
    r = length(bm.a)
    λs = append!(collect(coefficients(bm.Λ)), zeros(I, r-degree(bm.Λ)-1))
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
