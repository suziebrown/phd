### Attempt to make an illustrative plot showing the problem of ancestral
#   degeneracy within particle Gibbs.
function smc_example(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function, delta::Float64, sigma::Float64)
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

function csmc_example(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function, immortal_positions::Array{Float64,1}, delta::Float64, sigma::Float64)

    # pre-allocate memory
    parents = Array{Int64, 2}(undef, T, N)
    positions = Array{Float64, 2}(undef, T+1, N)
    weights = Array{Float64, 2}(undef, T+1, N)

    # initialise
    positions[1, :] = initialsam(N, delta, sigma)
    immortal_indices = rand(1:N, T+1) # pre-sample all of them to speed up computation
    positions[1, immortal_indices[1]] = immortal_positions[1]
    weights[1,:] = potential(positions[1,:], observations[1], delta, sigma)
    weights[1,:] = weights[1,:] / sum(weights[1,:]) # normalise weights

    for t in 1:T # note: index t+1 corresponds to generation t
        # choose parents
        parents[t, :] = sample(1:N, Weights(weights[t,:]), N)
        parents[t, immortal_indices[t+1]] = immortal_indices[t]

        # update positions
        positions[t+1,:] = outransition(positions[t, parents[t, :]], delta, sigma)
        positions[t+1, immortal_indices[t+1]] = immortal_positions[t+1]

        # compute weights
        weights[t+1, :] = potential(positions[t+1, :], observations[t+1], delta, sigma)
        weights[t+1, :] = weights[t+1,:] / sum(weights[t+1,:])
    end

    return (positions = positions, weights = weights, parents=parents, immortal=immortal_indices)
end

function trace_lineage(N::Int64, T::Int64, parents::Array{Int64,2})
    ## pre-allocate memory
    lineages = Array{Int64, 2}(undef, T+1, N)
    ## initialise
    lineages[T+1,:] = collect(1:N)

    for t in T:-1:1
        for i in 1:N
            lineages[t,i] = parents[t,lineages[t+1,i]]
        end
    end

    return lineages
end

## Generate data
delta = 0.05
sigma = 1.0
T=40
N=20
myobs = ousim(T, delta, sigma, false)

mypreviter = smc_example(N,T,myobs,ouinit,outransition,oupotential, delta, sigma)
myprevlineages = trace_lineage(N,T,mypreviter.parents)
imm_index = sample(1:N, Weights(mypreviter.weights[T,:]))
myimmindex = myprevlineages[:,imm_index]
myimmpos = Array{Float64,1}(undef,T+1)
for t in 1:(T+1)
    myimmpos[t] = mypreviter.positions[t,myimmindex[t]]
end
mysam = csmc_example(N,T,myobs,ouinit,outransition,oupotential,myimmpos, delta, sigma)


## Compute genealogy
gene_index = trace_lineage(N,T,mysam.parents)
gene_pos = Array{Float64,2}(undef,T+1,N)
for i in 1:N
    for t in 1:(T+1)
        gene_pos[t,i] = mysam.positions[t, gene_index[t,i]]
    end
end

## Plot genealogy
myancplot = plot(1:T+1, gene_pos[:,1], line=(:black), legend=false, xaxis="t", yaxis="position")#, ylims=(-0.5,N+0.5))
for i in 2:N
    plot!(1:T+1, gene_pos[:,i], line=:black)
end

# STILL NEED TO HIGHLIGHT THE IMMORTAL LINEAGE AND THE NEWLY SAMPLED LINEAGE

#title!("foo") ## change plot title for different cases
myancplot
#savefig("mylovelyplot.pdf") ## change file name for different cases
