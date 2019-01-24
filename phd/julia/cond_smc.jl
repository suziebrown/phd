# requires Random, Distributions, StatsBase

function csmc(N::Int64, T::Int64, observations::Array{Float64,1} initialsam::Function, transition::Function, potential::Function, immortal_positions::Array{Float64,1})
    # pre-allocate memory
    immortal_indices = Array{Int64, 1}(undef, T+1) # probably doesn't need initialising
    positions = Array{Float64, 1}(undef, N)  # probably doesn't need initialising
    # (also need to store ancestry tree but haven't worked out an efficient way to update it yet.)
    weights = Array{Float64,1}(undef, N)

    # set OU process parameters
    # should find a more flexible way to write this function, really for any model...
    delta = 1.0
    sigma = 0.1

    # initialise
    positions = initialsam(N, delta, sigma)
    immortal_indices = rand(1:N, T+1) # pre-sample all of them to speed up computation
    positions[immortal_indices[1]] = immortal_positions[1]
    for i in 1:N
        weights[i] = potential(positions[i], observations[1], delta, sigma) # VECTORISE ME
    end
    weights = weights./sum(weights) # normalise weights

    for t in 0:T-1 # indexing could be wrong because vectors start at 1...
        # choose parents
        parents = sample(1:N, Weights(weights), N)
        parents[immortal_indices[t+1]] = immortal_indices[t]
        # and store to tree structure

        # update positions
        positions = outransition(positions[parents], delta, sigma)
        positions[immortal_indices[t+1]] = immortal_positions[t+1]

        # compute weights
        for i in 1:N
            weights[i] = potential(positions[i], observations[t+1], delta, sigma) # VECTORISE ME
        end
        weights = weights./sum(weights) # normalise weights
    end
end
