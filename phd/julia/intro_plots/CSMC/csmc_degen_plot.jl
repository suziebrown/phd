### Make illustrative plot showing the problem of ancestral
#   degeneracy within particle Gibbs.
#
#   All of the required functions are contained in this file (I think!)
#

using Statistics, StatsBase, Random, Distributions, Plots

# generate a sequence of observations from an Orstein-Uhlenbeck process
function ousim(T::Int64, delta::Float64, sigma::Float64, returnstates::Bool)
    # pre-allocate arrays
    states = Array{Float64}(undef, T+1)
    observations = Array{Float64}(undef, T+1)

    # generate Normal noise:
    states_noise = rand(Normal(0, delta^(0.5)), T+1)
    obs_noise = rand(Normal(0, sigma), T+1)

    # initial states
    states[1] = states_noise[1]
    observations[1] = states[1] + obs_noise[1]

    # simulate the process
    for t in 2:T+1
        states[t] = (1-delta) * states[t-1] + states_noise[t]
        observations[t] = states[t] + obs_noise[t]
    end

    if returnstates
        return (states = states, observations = observations)
    else
        return observations
    end
end

# sample initial particle positions
function ouinit(N::Int64, delta::Float64, sigma::Float64)
    rand(Normal(0, delta^(0.5)), N)
end

# propagate states: sample new positions from current ones
function outransition(oldpos::Array{Float64,1}, delta::Float64, sigma::Float64)
    N = length(oldpos)
    rand(Normal(), N) .* delta^(0.5) .+ (1-delta) .* oldpos
end

# calculate potentials between particle positions and observations (to compute weights)
function oupotential(pos::Array{Float64,1}, obs::Float64, delta::Float64, sigma::Float64)
    (2*pi)^(-0.5) * sigma^(-1) * exp.(-(obs .- pos).^2 / (2*sigma^2))
end

# one iteration of SMC (for drawing the immortal trajectory)
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

# one iteration of CSMC (for drawing all the new trajectories)
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

# find the sequence of parental indices for each lineage
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

mypreviter = smc_example(N,T,myobs,ouinit,outransition,oupotential, delta, sigma) # run SMC to get ensemble of trajectories
myprevlineages = trace_lineage(N,T,mypreviter.parents) # track the ancetsral indices for each trajectory (actually only need it for the chosen trajectory)
imm_index = sample(1:N, Weights(mypreviter.weights[T,:])) # choose a lineage to keep as `prev iter' / 'immortal' output
myimmindex = myprevlineages[:,imm_index] # ancestral indices of chosen immortal lineage

# calculate positions of chosen immortal lineage
myimmpos = Array{Float64,1}(undef,T+1)
for t in 1:(T+1)
    myimmpos[t] = mypreviter.positions[t,myimmindex[t]]
end

# run csmc conditional on the chosen immortal lineage
mysam = csmc_example(N,T,myobs,ouinit,outransition,oupotential,myimmpos, delta, sigma)


## Compute genealogy positions
gene_index = trace_lineage(N,T,mysam.parents)
gene_pos = Array{Float64,2}(undef,T+1,N)
for i in 1:N
    for t in 1:(T+1)
        gene_pos[t,i] = mysam.positions[t, gene_index[t,i]]
    end
end

## Plot genealogy
myancplot = plot(1:T+1, gene_pos[:,1], line=(:black), legend=false, xaxis="t", yaxis="position", ticks = nothing, border = :none)#, ylims=(-0.5,N+0.5))
for i in 2:N
    plot!(1:T+1, gene_pos[:,i], line=:black) # add line for each trajectory
end
plot!(1:T+1, myimmpos, line=(5,:black)) # thick line shows immortal trajectory
nextiter_index = sample(1:N, Weights(mysam.weights[T,:])) # sample a new trajectory from csmc output
plot!(1:T+1, gene_pos[:,nextiter_index], line=(3,RGB(0.53, 0.0, 0.69))) # purple line shows newly sampled trajectory

myancplot # show plot
#savefig("PG_degen.pdf") # save plot
