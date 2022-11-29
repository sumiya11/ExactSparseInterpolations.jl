using BenchmarkTools

function fib(n)
    if n < 2
        return 1
    end
    return fib(n - 1) + fib(n - 2)
end

function gcdx1(a::Integer, b::Integer)
    T = promote_type(typeof(a), typeof(b))        
    # a0, b0 = a, b
    s0, s1 = oneunit(T), zero(T)
    t0, t1 = s1, s0
    # The loop invariant is: s0*a0 + t0*b0 == a   
    x = a % T
    y = b % T
    while y != 0
        q = div(x, y)
        x, y = y, rem(x, y)
        s0, s1 = s1, s0 - q*s1
        t0, t1 = t1, t0 - q*t1
    end
    x < 0 ? (-x, -s0, -t0) : (x, s0, t0)
end

function gcdx3(a::Integer, b::Integer)
    if b > a
        a, b = b, a
    end 
    T = promote_type(typeof(a), typeof(b))        
    # a0, b0 = a, b
    s0, s1 = oneunit(T), zero(T)
    t0, t1 = s1, s0
    # The loop invariant is: s0*a0 + t0*b0 == a   
    x = a % T
    y = b % T
    while y != 0
        q = div(x, y)
        x, y = y, x - q*y
        s0, s1 = s1, s0 - q*s1
        t0, t1 = t1, t0 - q*t1
    end
    x < 0 ? (-x, -s0, -t0) : (x, s0, t0)
end

function gcdx2(a::Integer, b::Integer)
    T = promote_type(typeof(a), typeof(b))        
    # a0, b0 = a, b
    s0, s1 = oneunit(T), zero(T)
    t0, t1 = s1, s0
    # The loop invariant is: s0*a0 + t0*b0 == a   
    x = a % T
    y = b % T
    while y != 0
        q = div(x, y, RoundNearestTiesAway)
        x, y = y, x - q*y
        s0, s1 = s1, s0 - q*s1
        t0, t1 = t1, t0 - q*t1
    end
    x < 0 ? (-x, -s0, -t0) : (x, s0, t0)
end


x, y = fib(29), fib(30)

@benchmark gcdx1($x, $y)

@benchmark gcdx2($x, $y)

@benchmark gcdx3($x, $y)
