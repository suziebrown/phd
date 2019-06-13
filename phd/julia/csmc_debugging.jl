function mrca_fullstore(ancestry::Array{Int64,2}, leaves::Array{Int64,1})
    T = size(ancestry)[1]
    mrca = T
    for t in 1:T
        leaves = ancestry[T-t+1, leaves]
        if allequal(leaves)
            mrca = t
            break
        end
    end
    return mrca
end

function smc_example(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function)
    # set OU process parameters
    delta = 1.0
    sigma = 0.1

    # pre-allocate memory
    parents = Array{Int64, 2}(undef, T, N)
    positions = Array{Float64, 2}(undef, T+1, N)
    weights = Array{Float64, 2}(undef, T+1, N)

    # initialise
    positions[1, :] = initialsam(N, delta, sigma)
    weights[1,:] = potential(positions[1,:], observations[1], delta, sigma)
    weights[1,:] = weights[1,:] / sum(weights[1,:]) # normalise weights

    for t in 1:T # note: index t+1 corresponds to generation t
        # choose parents
        parents[t, :] = sample(1:N, Weights(weights[t,:]), N)
        # update positions
        positions[t+1,:] = outransition(positions[t, parents[t, :]], delta, sigma)
        # compute weights
        weights[t+1, :] = potential(positions[t+1, :], observations[t+1], delta, sigma)
        weights[t+1, :] = weights[t+1,:] / sum(weights[t+1,:])
    end

    return (positions = positions, weights = weights, parents=parents)
end

N=10
T=10
testobs = ousim(T, 0.1, 0.1, false)
plot(testobs)
sam = smc_example(N,T,testobs,ouinit,outransition,oupotential)
println(sam.parents)

mrca_fullstore(sam.parents, collect(1:N))
