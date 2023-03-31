using SIMD
import SIMD: LVec
using BenchmarkTools
using Random
using LoopVectorization
import VectorizationBase
using StaticArrays
using Plots
using Test

struct StaticVect{A, B, C}

end

struct StaticVect2{V}

end

@generated function dotproduct(
        sv::StaticVect2{V}, 
        a::Vector{T}) where {V, T}
    s = :(0)
    for i in V.types
        if !iszero(i)
            s = :($s + $i * a[1])
        end
    end
    s
end

################################

"""
    merge_sorted_scalar!

    Merges two sorted arrays of integers `exps1` and `exps2` 
    into a single sorted array `buf_exps`. 

    Assumes `buf` has enough space to hold the result.

"""
function merge_sorted_scalar!(
        buf::Vector{E}, 
        exps1::Vector{E}, 
        exps2::Vector{E}) where {E}
    n1, n2 = length(exps1), length(exps2)
    i = 1
    i1, i2 = 1, 1
    @inbounds while i1 <= n1 && i2 <= n2
        if exps1[i1] < exps2[i2]
            buf[i] = exps1[i1]
            i1 = i1 + 1
        else
            buf[i] = exps2[i2]
            i2 = i2 + 1
        end
        i = i + 1
    end
    # Forgetting about the tails of the arrays for now
    # @inbounds while i1 <= n1
    #     buf_exps[i] = exps1[i1]
    #     i1 = i1 + 1
    #     i = i + 1
    # end
    # @inbounds while i2 <= n2
    #     buf_exps[i] = exps2[i2]
    #     i2 = i2 + 1
    #     i = i + 1
    # end
    buf
end

function merge_sorted_branchless!(
        buf::Vector{E},
        exps1::Vector{E},
        exps2::Vector{E}) where {E}
    n1, n2 = length(exps1), length(exps2)
    i = 1
    i1, i2 = 1, 1
    shift = E(sizeof(E) * 8 - 1)
    oneE = one(E)
    @inbounds while i1 <= n1 && i2 <= n2
        # flag = (exps1[i1] - exps2[i2]) >> shift
        flag = E(exps1[i1] < exps2[i2])
        buf[i] = flag*exps1[i1] + (oneE - flag)*exps2[i2]
        i1 = i1 + flag
        i2 = i2 + (1 - flag)
        i = i + 1
    end
    buf
end

function merge_sorted_heuristic!(
        buf::Vector{E}, 
        exps1::Vector{E}, 
        exps2::Vector{E}) where {E}
    n1, n2 = length(exps1) - 8, length(exps2) - 8
    i = 1
    i1, i2 = 1, 1
    @inbounds while i1 <= n1 && i2 <= n2
        if exps1[i1] < exps2[i2]
            if exps1[i1 + 7] < exps2[i2]
                mm1 = vload(Vec{8, E}, exps1, i1)
                vstore(mm1, buf, i)
                i1 = i1 + 8
                i = i + 8
            else
                buf[i] = exps1[i1]
                i1 = i1 + 1
                i = i + 1
            end
        else
            if exps1[i1] >= exps2[i2 + 7]
                mm1 = vload(Vec{8, E}, exps2, i2)
                vstore(mm1, buf, i)
                i2 = i2 + 8
                i = i + 8
            else
                buf[i] = exps2[i2]
                i2 = i2 + 1
                i = i + 1
            end
        end
    end
    buf
end

################################

function merge_sorted_simd_1!(
        buf::Vector{E},
        exps1::Vector{E},
        exps2::Vector{E}) where {E}
    n1, n2 = length(exps1) - 8, length(exps2) - 8
    i = 1
    i1, i2 = 1, 1
    @inbounds while i1 <= n1 && i2 <= n2
        a = vload(Vec{8, E}, exps1, i1)
        b = vload(Vec{8, E}, exps2, i2)
        _cmp = a - b
        _cmp = _cmp >> (sizeof(E) * 8 - 1)
        println("a = $a\nb = $b\n_cmp = $_cmp")
        k = -sum(_cmp)
        x = vblend(a, b, _cmp)
        vstore(x, buf, i)
        i1 += k
        i2 += 8 - k
        println("k = $k\nx = $x")
        println("i1 = $i1\ni2 = $i2\n")
        println("==========")
        i += 8
    end
    buf
end

################################

@generated function vmin(
        x::LVec{N, T}, 
        y::LVec{N, T},
        ) where {N, T <: SIMD.IntegerTypes}
    d = SIMD.Intrinsics.d
    llvmsuffix = SIMD.Intrinsics.suffix(N, T)
    s = """
    declare <$N x $(d[T])> @llvm.smin.$(llvmsuffix)(<$N x $(d[T])> %0, <$N x $(d[T])> %1)

    define <$N x $(d[T])> @entry(<$N x $(d[T])>, <$N x $(d[T])>) #0 {
    top:
        %res = call <$N x $(d[T])> @llvm.smin.$(llvmsuffix)(<$N x $(d[T])> %0, <$N x $(d[T])> %1)
        ret <$N x $(d[T])> %res
    }

    attributes #0 = { alwaysinline }
    """
    return :(
        $(Expr(:meta, :inline));
        Base.llvmcall(($s, "entry"), LVec{N, T}, Tuple{LVec{N, T}, LVec{N, T}}, x, y)
    )
end

@generated function vblend(
        x::LVec{N, T}, 
        y::LVec{N, T},
        mask::LVec{N, T}
        ) where {N, T <: SIMD.IntegerTypes}
    d = SIMD.Intrinsics.d
    s = """
    declare <32 x i8> @llvm.x86.avx2.pblendvb(<32 x i8> %0, <32 x i8> %1, <32 x i8> %2)

    define <$N x $(d[T])> @entry(<$N x $(d[T])> %0, <$N x $(d[T])> %1, <$N x $(d[T])> %2) #0 {
    top:
        %3 = bitcast <$N x $(d[T])> %0 to <32 x i8>
        %4 = bitcast <$N x $(d[T])> %1 to <32 x i8>
        %5 = bitcast <$N x $(d[T])> %2 to <32 x i8>
        %6 = call <32 x i8> @llvm.x86.avx2.pblendvb(<32 x i8> %3, <32 x i8> %4, <32 x i8> %5)
        %res = bitcast <32 x i8> %6 to <$N x $(d[T])> 
        ret <$N x $(d[T])> %res
    }

    attributes #0 = { alwaysinline }
    """
    return :(
        $(Expr(:meta, :inline));
        Base.llvmcall(($s, "entry"), LVec{N, T}, 
        Tuple{LVec{N, T}, LVec{N, T}, LVec{N, T}}, 
        x, y, mask
        )
    )
end

vmin(x::Vec{N, T}, y::Vec{N, T}) where {N, T} = Vec{N, T}(vmin(x.data, y.data))
vblend(x::Vec{N, T}, y::Vec{N, T}, mask::Vec{N, T}) where {N, T} = Vec{N, T}(vblend(x.data, y.data, mask.data))

a = SIMD.Vec{8, Int32}(19)
b = SIMD.Vec{8, Int32}((8, 7, 6, 0, 0, 2, 1, 0))
mask = SIMD.Vec{8, Int32}((0, 0, -1, -1, 0, -1, 0, -1))

mask = mask >> 31

vmin(a, b)
vblend(a, b, mask)

################################

function generate_data(E, n, chunksize=6:10)
    n1, n2 = n, n
    exps = unique(rand(E, n1 + n2))
    indices = rand(1:length(exps), div(length(exps), 2))
    exps1 = sort(exps[indices])
    exps2 = sort(exps[setdiff(1:length(exps), indices)[1:length(indices)]])
    n1, n2 = length(exps1), length(exps2)
    @assert n1 == n2
    lo, hi = 1, rand(chunksize)
    while hi <= min(n1, n2)
        coinflip = rand(Bool)
        chunk = sort(vcat(exps1[lo:hi], exps2[lo:hi]))
        k = div(length(chunk), 2)
        chunkmin = chunk[1:k]
        chunkmax = chunk[k+1:end]
        if coinflip
            exps1[lo:hi] = chunkmin
            exps2[lo:hi] = chunkmax
        else
            exps1[lo:hi] = chunkmax
            exps2[lo:hi] = chunkmin
        end  
        lo = hi + 1
        hi = hi + rand(chunksize)
    end
    buf = Vector{E}(undef, n1 + n2)
    exps1, exps2, buf
end

################################

E = Int32
n = 10
exps1, exps2, buf = generate_data(E, n, 20:30)
exps1, exps2 = sort(rand(E(1):E(2^20), n)), sort(rand(E(1):E(2^20), n))
buf = Vector{E}(undef, n + n)

A = merge_sorted_scalar!(deepcopy(buf), exps1, exps2)
B = merge_sorted_simd_1!(deepcopy(buf), exps1, exps2)


@btime merge_sorted_scalar!($buf, $exps1, $exps2);
@btime merge_sorted_branchless!($buf, $exps1, $exps2);
@btime merge_sorted_heuristic!($buf, $exps1, $exps2);
@btime merge_sorted_simd_1!($buf, $exps1, $exps2);

times_scalar = []
times_branchless = []
times_heuristic = []
begin
    E = UInt32
    n = 10^4
    chunksizes = [1:1, 10:20]
    for k in chunksizes
        chunksize = k
        exps1, exps2, buf = generate_data(UInt32, n, chunksize)
        
        stat = @benchmarkable merge_sorted_scalar!(buf, exps1, exps2)
        push!(times_scalar, minimum(run(stat, seconds=1.0)).time)

        stat = @benchmarkable merge_sorted_branchless!(buf, exps1, exps2)
        push!(times_branchless, minimum(run(stat, seconds=1.0)).time)

        stat = @benchmarkable merge_sorted_heuristic!(buf, exps1, exps2)
        push!(times_heuristic, minimum(run(stat, seconds=1.0)).time)
    end
end
begin
    p = plot(map(first, chunksizes), times_scalar, label="scalar", 
            title="$n elements of type $E",
            xaxis="Chunk size",
            yaxis="Time (ns)"
        )
    plot!(map(first, chunksizes), times_branchless, label="branchless")
    plot!(map(first, chunksizes), times_heuristic, label="heuristic")
end

# make batches of 100.
# graph dependent on batch size 1 -> 100 // 1000.

################################
################################

struct StaticVectorI{T}

end

function dotprod(
        sv::StaticVectorI{T}, 
        x) where {T}
    s = :()
    for (Y, i) in zip(T.types, 1:length(x))
        if !iszero(Y)
            if s == :()
                s = :($Y*x[$i])
            else
                s = :($s + $Y*x[$i])
            end
        end
    end
    :(@inbounds $s)
end

struct StaticVector4{A, B, C, D}

end

@generated function dotprod(
        sv::StaticVector4{A, B, C, D}, 
        x) where {A, B, C, D}
    s = :()
    for (Y, i) in zip((A, B, C, D), 1:4)
        if !iszero(Y)
            if s == :()
                s = :($Y*x[$i])
            else
                s = :($s + $Y*x[$i])
            end
        end
    end
    :(@inbounds $s)
end

################################
################################

function cmp_1(s, a, b)
    @inbounds for i in 1:length(s)
        if a[i] < b[i]
            s[i] = a[i]
        else
            s[i] = b[i]
        end
    end
    s
end
function cmp_2(s, a, b)
    @turbo for i in 1:length(s)
        s[i] = min(a[i], b[i])
    end
    s
end

n = 1000
a = rand(Int32, n)
b = rand(Int32, n)
s = similar(a)

@btime cmp_1($s, $a, $b);
@btime cmp_2($s, $a, $b);


