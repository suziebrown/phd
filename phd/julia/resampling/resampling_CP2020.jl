### RESAMPLING ALGORITHMS a la Chopin & Papaspiliopoulos
#   (based on ALgms 9.1--9.6 in their 2020 book)
using Distributions

function icdf(N::Int64, weight::Array{Float64,1}, seeds::Array{Float64,1})
    # seeds must be an ORDERED vector length N of numbers in [0,1]
    # weights must be a vector length N of numbers in [0,1] summing to 1
    # produces ordered indices (does not satisfy standing assumption!)

    # allocate memory:
    parents = Array{Float64,1}(undef, N)
    # initialise:
    s = weight[1]
    m = 1

    for i in 1:N
        while s < seeds[i]
            m = m + 1
            s = s + weight[m]
        end
        parents[i] = m
    end
    return parents
end

function seedmn(N::Int64)
    E = rand(Exponential(1), N+1)
    S = cumsum(E)
    S = S / S[N+1]
    return S[1:N]
end

function seedstrat(N::Int64)
    U = Array{Float64,1}(undef,N)
    for i in 1:N
        U[i] = rand(Uniform((i-1)/N, i/N))
    end
    return U
end

function seedsyst(N::Int64)
    # it's possible to implement a faster inverse-CDF method based on syst structure (see Alg 9.6 in book)
    U = Array{Float64,1}(undef,N)
    u = rand(Uniform(0,1/N))
    for i in 1:N
        U[i] = u + (i-1)/N
    end
    return U
end


## altogether now...

function resam_mn(N::Int64, weight::Array{Float64,1})
    Int.(icdf(N, weight, seedmn(N)))
end
function resam_strat(N::Int64, weight::Array{Float64,1})
    Int.(icdf(N, weight, seedstrat(N)))
end
function resam_syst(N::Int64, weight::Array{Float64,1})
    Int.(icdf(N, weight, seedsyst(N)))
end
function resam_star(N::Int64, weight::Array{Float64,1})
    parent = rand(Categorical(weight))
    parents = fill(parent, N)
    return (parents)
end
function resam_resmn(N::Int64, weight::Array{Float64,1})
    parents = Array{Int64,1}(undef,N)
    j=1
    floors = Int.(floor.(N.*weight))
    resids = (N.*weight - floors)
    R = Int(round(sum(resids)))
    resids = resids/R
    for i in 1:N # deterministic part
        parents[j:(j+floors[i]-1)] = repeat([i], floors[i])
        j = j + floors[i]
    end
    parents[(N-R+1):N] = icdf(R, resids, seedmn(R))
    return sort(parents)
end
function resam_resstrat(N::Int64, weight::Array{Float64,1})
    parents = Array{Int64,1}(undef,N)
    j=1
    floors = Int.(floor.(N.*weight))
    resids = (N.*weight - floors)
    R = Int(round(sum(resids)))
    resids = resids/R
    for i in 1:N # deterministic part
        parents[j:(j+floors[i]-1)] = repeat([i], floors[i])
        j = j + floors[i]
    end
    parents[(N-R+1):N] = icdf(R, resids, seedstrat(R))
    return sort(parents)
end
function resam_ressyst(N::Int64, weight::Array{Float64,1})
    parents = Array{Int64,1}(undef,N)
    j=1
    floors = Int.(floor.(N.*weight))
    resids = (N.*weight - floors)
    R = Int(round(sum(resids)))
    resids = resids/R
    for i in 1:N # deterministic part
        parents[j:(j+floors[i]-1)] = repeat([i], floors[i])
        j = j + floors[i]
    end
    parents[(N-R+1):N] = icdf(R, resids, seedsyst(R))
    return sort(parents)
end
function resam_resstar(N::Int64, weight::Array{Float64,1})
    parents = Array{Int64,1}(undef,N)
    j=1
    floors = Int.(floor.(N.*weight))
    resids = (N.*weight - floors)
    R = Int(round(sum(resids)))
    resids = resids/R
    for i in 1:N # deterministic part
        parents[j:(j+floors[i]-1)] = repeat([i], floors[i])
        j = j + floors[i]
    end
    resparent = rand(Categorical(resids))
    parents[(N-R+1):N] = fill(resparent, R)
    return Int.(sort(parents))
end
function resam_ssp(N::Int64, weight::Array{Float64,1})
    parents = N.*weight
    n = 1
    m = 2

    while n<N || m<N
        println(n)
        println(m)
        println(parents[n])
        println(parents[m])

        delta_n = floor(parents[n])+1 - parents[n]
        delta_m = parents[m] - floor(parents[m])
        delta = min(delta_n, delta_m)
        eps_n = parents[n] - floor(parents[n])
        eps_m = floor(parents[m])+1 - parents[m]
        eps = min(eps_n, eps_m)

        u = rand()
        if u <= eps/(eps+delta) #  then use delta
            if delta_n <= delta_m # then delta==delta_n
                parents[n] = floor(parents[n]) +1
                parents[m] = parents[m] - delta
                n = max(n,m)+1
            else # then delta==delta_m
                parents[n] = parents[n] + delta
                parents[m] = floor(parents[m])
                m = max(n,m)+1
            end
        else # then use eps
            if eps_n <= eps_m # then eps==eps_n
                parents[n] = floor(parents[n])
                parents[m] = parents[m] + eps
                n = max(n,m)+1
            else # then eps==eps_m
                parents[n] = parents[n] - eps
                parents[m] = floor(parents[m]) +1
                m = max(n,m)+1
            end
        end

        println(parents[n])
        println(parents[m])
        println()
    end

    return Int.(parents)
end
function resam_ssp2(N::Int64, weight::Array{Float64,1})
    parents = N.*weight
    n = 1
    m = 2
    iters=0

    while n<N || m<N
        println(n)
        println(m)
        println(parents[n])
        println(parents[m])

        delta = min(floor(parents[n])+1 - parents[n], parents[m] - floor(parents[m]))
        eps = min(parents[n] - floor(parents[n]), floor(parents[m])+1 - parents[m])

        u = rand()
        if u <= eps/(eps+delta) #  then use delta
            parents[m] = parents[m] - delta
            parents[n] = parents[n] + delta
            if parents[m]==floor(parents[m])
                if parents[n]==floor(parents[n])
                    n=m+1
                    m=m+2
                else
                    m=m+1
                end
            elseif parents[n]==floor(parents[n])
                n=m
                m=m+1
            end
        else # then use eps
            parents[m] = parents[m] + eps
            parents[n] = parents[n] - eps
            if parents[m]==floor(parents[m])
                if parents[n]==floor(parents[n])
                    n=m+1
                    m=m+2
                else
                    m=m+1
                end
            elseif parents[n]==floor(parents[n])
                n=m+1
            end
        end

        iters = iters+1
        if iters>100
            break
        end

        println(parents[n])
        println(parents[m])
        println()
    end

    return parents
end
## TEST
N=20
weight = rand(Exponential(),N)
weight = weight / sum(weight)
# println(round.(weight, digits=4))
# println( resam_mn(N,weight) )
# println( resam_strat(N,weight) )
# println( resam_syst(N,weight) )
# println( resam_star(N,weight) )
# println( resam_resmn(N,weight) )
# println( resam_resstrat(N,weight) )
# println( resam_ressyst(N,weight) )
# println( resam_resstar(N,weight) )
println( resam_ssp2(N,weight) )
