### RESAMPLING ALGORITHMS a la Chopin & Papaspiliopoulos
#   (based on ALgms 9.1--9.6 in their 2020 book)

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

## TEST
# N=20
# weight = rand(Exponential(),N)
# weight = weight / sum(weight)
# println(round.(weight, digits=4))
# println( resam_mn(N,weight) )
# println( resam_strat(N,weight) )
# println( resam_syst(N,weight) )
# println( resam_resmn(N,weight) )
# println( resam_resstrat(N,weight) )
# println( resam_ressyst(N,weight) )
