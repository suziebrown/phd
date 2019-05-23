# requires Random, Distributions, StatsBase

function csmc_fullstore(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function, immortal_positions::Array{Float64,1})
    # set OU process parameters
    # (should find a more flexible way to write this function, really for any model...)
    # (allow a (named) vector (any length) of 'parameters' as input)
    delta = 1.0
    sigma = 0.1

    # INITIALISE ANCESTRAL TREE - but how??
    parents = Array{Int64, 2}(undef, T, N)

    # initialise
    positions = initialsam(N, delta, sigma)
    immortal_indices = rand(1:N, T+1) # pre-sample all of them to speed up computation
    positions[immortal_indices[1]] = immortal_positions[1]
    weights = potential(positions, observations[1], delta, sigma)
    weights = weights / sum(weights) # normalise weights

    for t in 1:T # note: index t+1 corresponds to generation t
        # choose parents
        parents[t, :] = sample(1:N, Weights(weights), N)
        parents[t, immortal_indices[t+1]] = immortal_indices[t]
        # STORE TO ANCESTRAL TREE

        # update positions
        positions = outransition(positions[parents[t, :]], delta, sigma)
        positions[immortal_indices[t+1]] = immortal_positions[t+1]

        # compute weights
        weights = potential(positions, observations[t+1], delta, sigma)
        weights = weights / sum(weights)
    end

    return (parents=parents, immortal=immortal_indices)
    # RETURN ANCESTRAL TREE
end


#= building ancestral tree...

Can use
    findall(x->x==i, parents)
to return the children indices of a particular parent i.
But idk how to add these as TreeNodes to parent.children in a vectorised way.
Maybe just use a loop for now.
Concatenate arrays using
    vcat(A,B)
:)
=#
