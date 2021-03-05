### PLOT A REALISATION OF tau_N(t) vs. t
#
# Project Objectives:
# 1. compare tau_N for different resampling schemes
# 2. explore behaviour of tau_N as T,N vary

using Random, Distributions, StatsBase, Plots, Statistics
# ALSO LOAD ALL FUNCTIONS FROM:
# phd/julia/OU_model.jl
# phd/julia/resampling/resampling_CP2020.jl

function smc_example(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function)
    # pre-allocate memory
    parents = Array{Int64, 2}(undef, T, N)
    cN = Array{Float64,1}(undef, T)
    positions = Array{Float64, 2}(undef, T+1, N)
    weights = Array{Float64, 2}(undef, T+1, N)
    offsprcounts = Array{Int64,1}(undef, N) # for temporary storage within loop

    # initialise
    positions[1, :] = initialsam(N, delta, sigma)
    weights[1,:] = potential(positions[1,:], observations[1], delta, sigma)
    weights[1,:] = weights[1,:] / sum(weights[1,:]) # normalise weights

    # pre-calculate random numbers for OU transitions
    rands = reshape(rand(Normal(), N*T), T, N)

    for t in 1:T # note: index t+1 corresponds to generation t
        # choose parents
        parents[t, :] = sample(1:N, Weights(weights[t,:]), N)
        # calculate offspring counts
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents[t, :])
        end
        # calculate cN
        cN[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        # update positions
        positions[t+1,:] = outransition(positions[t, parents[t, :]], delta, sigma, rands[t, :])
        # compute weights
        weights[t+1, :] = potential(positions[t+1, :], observations[t+1], delta, sigma)
        weights[t+1, :] = weights[t+1,:] / sum(weights[t+1,:])
    end

    return cN
end

### Let's go!

T=10 #time horizon
N=20 # number of particles
# OU process parameters:
delta = 1.0 # step size & transition noise (variance)
sigma = 0.1 # observation noise (s.d.)
# generate realisation
mysimobs = ousim(T, delta, sigma, false)
mysimcN = smc_example(N, T, mysimobs, ouinit, outransition, oupotential)
sumcN = cumsum(mysimcN)

# PLOTS
# plot observations:
scatter(0:T, mysimobs, marker=(:circle, 0.6, :black), leg=false)
# plot tau_N:
plot([0;sumcN], 0:T, line=(:steppre, 0.6, :black), leg=false, xaxis="t", yaxis="tau_N(t)")
