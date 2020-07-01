function S61_SIEVE(n::Integer)
    done = Threads.Atomic{Bool}(false)
    primes = Int[]
    sieves = [Channel{Int}() for _ = 1:n]
    for i in 1:n
        Threads.@spawn begin
            mp = p = take!(sieves[i])
            push!(primes, p)
            if i == n
                done[] = true
                return
            end
            for m in sieves[i]
                while m > mp; mp += p; end
                if m < mp
                    put!(sieves[i + 1], m)
                end
            end
        end
    end

    put!(sieves[1], 2)
    k = 3
    while !done[]
        put!(sieves[1], k)
        k += 2
    end
    return primes
end

primes = S61_SIEVE(100)
println("Found $(length(primes)) primes:")
println(primes)
