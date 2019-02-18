## make motivational plots
using Random, Distributions, StatsBase, Plots, Statistics
# small values of N,T only!
function csmc_example(N::UInt16, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function, immortal_positions::Array{Float64,1})
    # set OU process parameters
    delta = 1.0
    sigma = 0.1

    # pre-allocate memory
    parents = Array{UInt16, 2}(undef, T, N)
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

function smc_example(N::UInt16, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function)
    # set OU process parameters
    delta = 1.0
    sigma = 0.1

    # pre-allocate memory
    parents = Array{UInt16, 2}(undef, T, N)
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

##--- required functions ---

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
function ouinit(N::UInt16, delta::Float64, sigma::Float64)
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


##--- run example ---

# set constants
delta = 0.1 # noise variance in AR(1) process
sigma = 0.1 # noise s.d. in observations
T = Int64(23) # number of generations/time steps
N = UInt16(20) # total number of particles

# generate observations & immortal trajectory
observations = ousim(T, delta, sigma, false)
rtsout = ourts(delta, sigma, observations) # .+ rand(Normal(0,0.1),T+1)
kalout = oukalman(delta,sigma,observations)
smcout = smc_example(N, T, observations, ouinit, outransition, oupotential)


##--- make plots ---
scatter(smcout.positions, color=:black, ms=smcout.weights*N/2, grid=:x, xaxis=("time", (0,T+2), 1:1:T+1), legend=false, axis=("position", (-1.25,0.5)), size=(600,250))
plot!(rtsout.mean, color=:purple, ribbon=(2*rtsout.variance .^0.5, 2*rtsout.variance .^0.5), fillalpha=0.3)


##--- save most recent plot ---
savefig("smc_kalman_2.pdf")
