using StatsBase, Random, Distributions, Plots
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

function oukalman(delta::Float64, sigmain::Float64, observations::Array{Float64,1})
    T = length(observations) - 1

    # assign memory
    xhat = Array{Float64, 1}(undef, T+1) # means \hat{x}_{k|k}
    sigma = Array{Float64,1}(undef, T+1) # variances \Sigma_{k|k}

    # initialise
    tempx = 0 # \hat{x}_{k+1|k}
    tempsigma = 1 # \Sigma_{k+1|k}

    # forward recursion
    for t in 1:T
        a = tempsigma / (tempsigma + sigmain^2)
        xhat[t] = tempx + a * (observations[t] - tempx)
        sigma[t] = tempsigma * (1 - a)
        tempx = (1-delta) * xhat[t]
        tempsigma = (1-delta)^2 * sigma[t] + delta
    end

    # final filter state
    a = tempsigma / (tempsigma + sigmain^2)
    xhat[T+1] = tempx + a * (observations[T+1] - tempx)
    sigma[T+1] = tempsigma * (1 - a)

    return (mean = xhat, variance = sigma)
end

# could be useful to have a different version that takes kalman output as input...
function ourts(delta::Float64, sigmain::Float64, observations::Array{Float64,1})
    # forward Kalman pass
    kalman = oukalman(delta, sigmain, observations)
    xhat = kalman.mean
    sigma = kalman.variance
    # Warning: this function will incrementally OVERWRITE this Kalman output
    # (xhat, sigma) with the RTS output

    # backward recursion
    for i in 1:T
        t = T+1 - i # reverse time
        a = (1-delta) * sigma[t] / ((1-delta)^2 * sigma[t] + delta)
        xhat[t] = xhat[t] + a * (xhat[t+1] - (1-delta) * xhat[t])
        sigma[t] = sigma[t] + a^2 * (sigma[t+1] - (1-delta) * sigma[t] - delta)
    end

    return (mean = xhat, variance = sigma)
end

allequal(x) = all(y->y==x[1],x)

# should be called treeheight or something really
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

#----------------
# HERE BEGINS THE SIMULATION
#----------------

# set constants
delta = 0.1 # noise variance in AR(1) process
sigma = 0.1 # noise s.d. in observations
T = Int64(10000) # number of generations/time steps
N = Int64(8192) # total number of particles

# generate observations & immortal trajectory
observations = ousim(T, delta, sigma, false)
immpos = ourts(delta, sigma, observations).mean

##--- subtree sampling

# set constants
nvals = [2,4,8,16,32,64,128,256,512,1024,2048,4096,8192] # number of leaves in sampled subtree
nrep = 100

# initialise local variables
height = Array{Int64, 1}(undef, nrep)

# initialise output variables
meanall = Array{Float64, 1}(undef, length(nvals))
sdall = Array{Float64, 1}(undef, length(nvals))
noob = zeros(Int64, length(nvals))

for j in 1:length(nvals)
    println("now starting reps for n=", nvals[j])
    # uniformly sample 'nrep'x subtrees of size 'nvals[j]' & calculate height
    for i in 1:nrep
        # report progress (only works for nrep being multiple of 10)
        if (i*10) % nrep ==0
            println(100*i/nrep, "% of reps completed")
        end

        # run csmc
        csmcout = csmc_fullstore(N, T, observations, ouinit, outransition, oupotential, immpos)
        # sample a subtree
        leaves = sample(1:N, nvals[j])
        height[i] = mrca_fullstore(csmcout.parents, leaves)

        # error catching
        noob[j] += sum(height==T)
    end
    # mean & SD of tree heights
    meanall[j] = mean(height)
    sdall[j] = std(height)
end

# catching error of T being too small for N, i.e. reports no. of cases where treeheight was T
println("number of cases hitting limit was ", sum(noob))
println("means were: ", meanall)
println("SDs were: ", sdall)
# plot output (ribbon shows +/- 1 standard error)
plot(nvals, meanall/N, ribbon=(sdall*nrep^(-0.5)/N, sdall*nrep^(-0.5)/N), fill=:purple, fillalpha=0.25, leg=false, xaxis=:log10, line=(:purple), marker=(:purple), markerstrokecolor=:purple, title="CSMC treeheight, immortal=MAP, N=$N", xlabel="n", ylabel="average tree height /N")
