
# Convenient blackbox representation for polynomials and functions
mutable struct Blackbox{F}
    # the function to be evaluated
    f::F
    # the history of evaluations: x points and y points
    xs::Vector{Any}
    ys::Vector{Any}
end

Blackbox(f) = Blackbox(f, Any[], Any[])
Blackbox(f::Nemo.RingElem) = Blackbox(x -> evaluate(f, x))

function Blackbox(pq::AbstractAlgebra.Generic.Frac{T}) where {T<:Nemo.RingElem}
    Blackbox(x -> evaluate(numerator(pq), x) // evaluate(denominator(pq), x))
end

function Blackbox(pq::Tuple{T, T}) where {T<:Nemo.RingElem}
    Blackbox(x -> evaluate(pq[1], x) // evaluate(pq[2], x))
end

xs(bb::Blackbox) = bb.xs
ys(bb::Blackbox) = bb.ys

Base.count(bb::Blackbox) = (@assert length(bb.ys) == length(bb.xs); length(bb.ys))
Base.empty!(bb::Blackbox) = (empty!(bb.xs); empty!(bb.ys); bb)

# evaluate the blackbox
function (bb::Blackbox)(x)
    y = bb.f(x)
    push!(bb.xs, x)
    push!(bb.ys, y)
    y
end
