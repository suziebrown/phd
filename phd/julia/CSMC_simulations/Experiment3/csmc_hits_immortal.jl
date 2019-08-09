using Statistics, StatsBase, Random, Distributions, Plots, Dates
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
function mrca_hitimm(ancestry::Array{Int64,2}, leaves::Array{Int64,1}, immortal::Array{Int64,1})
    T = size(ancestry)[1]
    mrca = T
    hitimm = false
    for t in 1:T
        leaves = ancestry[T-t+1, leaves]
        if hitimm==false
            if immortal[T-t+1] in leaves
                hitimm=true
            end
        end
        if allequal(leaves)
            mrca = t
            break
        end
    end
    return (mrca=mrca, hitimm=hitimm)
end

#----------------
# HERE BEGINS THE SIMULATION
#----------------

# set constants
delta = 0.1 # noise variance in AR(1) process
sigma = 0.1 # noise s.d. in observations
N = Int64(1024) # total number of particles
T = Int64(150*N) # number of generations/time steps
nvals = [2,4,8,16,32,64,128,256,512,1023] # values of n to sample with (n=N not allowed)
nrep = 100 # number of repetitions for each n value

# generate observations & immortal trajectory
observations = ousim(T, delta, sigma, false)
imm = ourts(delta, sigma, observations)

# set less-constants
nsd = 0
immpos = imm.mean + nsd*(imm.variance).^(0.5)

# initialise local variables
heighttrue = Array{Int64, 1}(undef, nrep)
hitimm = Array{Bool, 1}(undef, nrep)
height1 = Array{Float64, 1}(undef, nrep)
height0 = Array{Float64, 1}(undef, nrep)
# initialise output variables
nsam0 = Array{Float64, 1}(undef, length(nvals))
nsam1 = Array{Float64, 1}(undef, length(nvals))
mean0 = Array{Union{Missing, Float64}, 1}(undef, length(nvals))
mean1 = Array{Union{Missing, Float64}, 1}(undef, length(nvals))
lquant0 = Array{Union{Missing, Float64}, 1}(undef, length(nvals))
uquant0 = Array{Union{Missing, Float64}, 1}(undef, length(nvals))
lquant1 = Array{Union{Missing, Float64}, 1}(undef, length(nvals))
uquant1 = Array{Union{Missing, Float64}, 1}(undef, length(nvals))
noob = zeros(Int64, length(nvals))

for j in 1:length(nvals)
    # report progress
    println("Starting n=", nvals[j], " at ", Dates.now())

    # uniformly sample 'nrep'x subtrees of size 'nvals[j]' & calculate height
    for i in 1:nrep

        # run csmc
        csmcout = csmc_fullstore(N, T, observations, ouinit, outransition, oupotential, immpos)
        # sample a subtree not containing immortal leaf
        immpart = csmcout.immortal[T+1]
        nonimm = deleteat!(collect(1:N), immpart)
        leaves = sample(nonimm, nvals[j], replace=false)
        # calculate tree heights
        mrcaout = mrca_hitimm(csmcout.parents, leaves, csmcout.immortal)
        heighttrue[i] = mrcaout.mrca
        hitimm[i] = mrcaout.hitimm

        # error catching
        noob[j] += (heighttrue[i] >= T)
        # report progress (only works for nrep being multiple of 10)
        if (i*10) % nrep ==0
            println(100*i/nrep, "% of reps completed")
        end
    end
    # normalise by N
    height0 = heighttrue[.!hitimm]/N
    height1 = heighttrue[hitimm]/N
    # statistics for tree heights (immortal line excluded):
    nsam0[j] = sum(.!hitimm)
    if nsam0[j]>0
        mean0[j] = mean(height0)
        lquant0[j] = quantile!(height0, 0.05)
        uquant0[j] = quantile!(height0, 0.95, sorted=true)
    else
        mean0[j] = missing
        lquant0[j] = missing
        uquant0[j] = missing
    end
    # statistics for tree heights (immortal line included):
    nsam1[j] = sum(hitimm)
    if nsam1[j]>0
        mean1[j] = mean(height1)
        lquant1[j] = quantile!(height1, 0.05)
        uquant1[j] = quantile!(height1, 0.95, sorted=true)
    else
        mean1[j] = missing
        lquant1[j] = missing
        uquant1[j] = missing
    end
end

# catching error of T being too small for N, i.e. reports no. of cases where treeheight was T
println("number of cases hitting limit was ", sum(noob))

# save results to file
datetime = Dates.now()
open("results_exp3_6", "w") do f
    write(f, "Simulation 6 for CSMC controlling for coalescing with immortal line \n
        File written at $datetime \n
        OU process with delta=$delta and sigma=$sigma \n
        Immortal line = MAP + $nsd SD \n
        T=$T, N=$N, reps per n =$nrep \n
        Values of n:\n
        $nvals \n
        Number of runs hitting limit:\n
        $noob \n\n
        ---Samples not hitting immortal line---\n
        Number of samples:\n
        $nsam0 \n
        Mean tree height/N:\n
        $mean0 \n
        Lower quantiles of tree height/N:\n
        $lquant0 \n
        Upper quantiles of tree height/N:\n
        $uquant0 \n
        ---Samples coalescing with immortal line---\n
        Number of samples:\n
        $nsam1 \n
        Mean tree height/N:\n
        $mean1 \n
        Lower quantiles of tree height/N:\n
        $lquant1 \n
        Upper quantiles of tree height/N:\n
        $uquant1 \n
        Sampled observation sequence:\n
        $observations \n"
    )
end

#---- plot results ----
plot()
plot!(nvals, mean0, xaxis=:log10, line=:dot, label="0", seriescolor=:purple)
plot!(nvals, mean1, label="1", seriescolor=:purple)

# proportion of samples coalescing with immortal line:
plot(nvals, nsam1/nrep, ylims=(0,1), xaxis=:log10, label="MAP")
plot!(nvals, nsam1/nrep, label="MAP")

#savefig("proportion_coalescing_immortal.pdf")

## NOte :: I accidentally generated new obs for nsd=1. I am using those going
## forward, but the nsd=0 results are made with different observations.
