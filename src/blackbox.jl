
mutable struct Blackbox{F}
    f::F
    xs::Vector{Any}
    ys::Vector{Any}
end

function Blackbox(f::Function)
    Blackbox(f, [], [])
end

function Blackbox(f::Nemo.RingElem)
    Blackbox(x -> evaluate(f, x), [], [])
end

function Blackbox(pq::AbstractAlgebra.Generic.Frac{T}) where {T<:Nemo.RingElem}
    Blackbox(x -> evaluate(numerator(pq), x) // evaluate(denominator(pq), x), [], [])
end

function Blackbox(pq::Tuple{T, T}) where {T<:Nemo.RingElem}
    Blackbox(x -> evaluate(pq[1], x) // evaluate(pq[2], x), [], [])
end

function xs(bb::Blackbox)
    bb.xs
end

function ys(bb::Blackbox)
    bb.ys
end

function Base.count(bb::Blackbox)
    length(bb.ys)
end

function Base.empty!(bb::Blackbox)
    empty!(bb.xs)
    empty!(bb.ys)
    bb
end

function (bb::Blackbox)(x)
    y = bb.f(x)
    push!(bb.xs, x)
    push!(bb.ys, y)
    y
end

