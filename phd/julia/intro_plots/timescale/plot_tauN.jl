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
        parents[t, :] = resam_mn(N, weights[t,:])
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

# Now a memory-efficient implementation that doesn't save stuff we don't care about
function smc_example_fast(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function)
    # pre-allocate memory
    parents = Array{Int64, 1}(undef, N)
    cN = Array{Float64,1}(undef, T) # for output
    positions = Array{Float64, 1}(undef, N)
    weights = Array{Float64, 1}(undef, N)
    offsprcounts = Array{Int64,1}(undef, N) # for temporary storage within loop

    # initialise
    positions = initialsam(N, delta, sigma)
    weights = potential(positions, observations[1], delta, sigma)
    weights = weights / sum(weights) # normalise weights

    # pre-calculate random numbers for OU transitions
    rands = reshape(rand(Normal(), N*T), T, N)

    for t in 1:T # note: index t+1 corresponds to generation t
        # choose parents
        parents = resam_mn(N, weights)
        # calculate offspring counts
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents)
        end
        # calculate cN
        cN[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        # update positions
        positions = outransition(positions[parents], delta, sigma, rands[t, :])
        # compute weights
        weights = potential(positions, observations[t+1], delta, sigma)
        weights = weights / sum(weights)
    end

    return cN
end

# ** Main function ** Now with various different resampling schemes, but the same transition rands
function smc_example_compare(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function)
    # pre-allocate memory
    parents_mn = Array{Int64, 1}(undef, N)
    cN_mn = Array{Float64,1}(undef, T) # for output
    positions_mn = Array{Float64, 1}(undef, N)
    weights_mn = Array{Float64, 1}(undef, N)
    parents_strat = Array{Int64, 1}(undef, N)
    cN_strat = Array{Float64,1}(undef, T) # for output
    positions_strat = Array{Float64, 1}(undef, N)
    weights_strat = Array{Float64, 1}(undef, N)
    parents_syst = Array{Int64, 1}(undef, N)
    cN_syst = Array{Float64,1}(undef, T) # for output
    positions_syst = Array{Float64, 1}(undef, N)
    weights_syst = Array{Float64, 1}(undef, N)
    parents_resmn = Array{Int64, 1}(undef, N)
    cN_resmn = Array{Float64,1}(undef, T) # for output
    positions_resmn = Array{Float64, 1}(undef, N)
    weights_resmn = Array{Float64, 1}(undef, N)
    parents_resstrat = Array{Int64, 1}(undef, N)
    cN_resstrat = Array{Float64,1}(undef, T) # for output
    positions_resstrat = Array{Float64, 1}(undef, N)
    weights_resstrat = Array{Float64, 1}(undef, N)
    parents_ressyst = Array{Int64, 1}(undef, N)
    cN_ressyst = Array{Float64,1}(undef, T) # for output
    positions_ressyst = Array{Float64, 1}(undef, N)
    weights_ressyst = Array{Float64, 1}(undef, N)
    offsprcounts = Array{Int64,1}(undef, N) # for temporary storage within loop

    # initialise
    positions_mn = initialsam(N, delta, sigma)
    weights_mn = potential(positions_mn, observations[1], delta, sigma)
    weights_mn = weights_mn / sum(weights_mn) # normalise weights
    positions_strat = initialsam(N, delta, sigma)
    weights_strat = potential(positions_strat, observations[1], delta, sigma)
    weights_strat = weights_strat / sum(weights_strat)
    positions_syst = initialsam(N, delta, sigma)
    weights_syst = potential(positions_syst, observations[1], delta, sigma)
    weights_syst = weights_syst / sum(weights_syst)
    positions_resmn = initialsam(N, delta, sigma)
    weights_resmn = potential(positions_resmn, observations[1], delta, sigma)
    weights_resmn = weights_resmn / sum(weights_resmn)
    positions_resstrat = initialsam(N, delta, sigma)
    weights_resstrat = potential(positions_resstrat, observations[1], delta, sigma)
    weights_resstrat = weights_resstrat / sum(weights_resstrat)
    positions_ressyst = initialsam(N, delta, sigma)
    weights_ressyst = potential(positions_ressyst, observations[1], delta, sigma)
    weights_ressyst = weights_ressyst / sum(weights_ressyst)

    # pre-calculate random numbers for OU transitions
    rands = reshape(rand(Normal(), N*T), T, N)

    for t in 1:T # note: index t+1 corresponds to generation t
        # choose parents
        parents_mn = resam_mn(N, weights_mn)
        # calculate offspring counts
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents_mn)
        end
        # calculate cN
        cN_mn[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        # update positions
        positions_mn = outransition(positions_mn[parents_mn], delta, sigma, rands[t, :])
        # compute weights
        weights_mn = potential(positions_mn, observations[t+1], delta, sigma)
        weights_mn = weights_mn / sum(weights_mn)

        parents_strat = resam_strat(N, weights_strat)
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents_strat)
        end
        cN_strat[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        positions_strat = outransition(positions_strat[parents_strat], delta, sigma, rands[t, :])
        weights_strat = potential(positions_strat, observations[t+1], delta, sigma)
        weights_strat = weights_strat / sum(weights_strat)

        parents_syst = resam_syst(N, weights_syst)
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents_syst)
        end
        cN_syst[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        positions_syst = outransition(positions_syst[parents_syst], delta, sigma, rands[t, :])
        weights_syst = potential(positions_syst, observations[t+1], delta, sigma)
        weights_syst = weights_syst / sum(weights_syst)

        parents_resmn = resam_resmn(N, weights_resmn)
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents_resmn)
        end
        cN_resmn[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        positions_resmn = outransition(positions_resmn[parents_resmn], delta, sigma, rands[t, :])
        weights_resmn = potential(positions_resmn, observations[t+1], delta, sigma)
        weights_resmn = weights_resmn / sum(weights_resmn)

        parents_resstrat = resam_resstrat(N, weights_resstrat)
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents_resstrat)
        end
        cN_resstrat[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        positions_resstrat = outransition(positions_resstrat[parents_resstrat], delta, sigma, rands[t, :])
        weights_resstrat = potential(positions_resstrat, observations[t+1], delta, sigma)
        weights_resstrat = weights_resstrat / sum(weights_resstrat)

        parents_ressyst = resam_ressyst(N, weights_ressyst)
        for i in 1:N
            offsprcounts[i] = count(j->(j==i), parents_ressyst)
        end
        cN_ressyst[t] = sum(offsprcounts .* (offsprcounts .- 1))./(N*(N-1))
        positions_ressyst = outransition(positions_ressyst[parents_ressyst], delta, sigma, rands[t, :])
        weights_ressyst = potential(positions_ressyst, observations[t+1], delta, sigma)
        weights_ressyst = weights_ressyst / sum(weights_ressyst)
    end

    return (cN_mn=cN_mn, cN_strat=cN_strat, cN_syst=cN_syst, cN_resmn=cN_resmn, cN_resstrat=cN_resstrat, cN_ressyst=cN_ressyst)
end

### Let's go!

T=150 #time horizon
N=150 # number of particles
# OU process parameters:
# sigma >> delta means observations are less informative, so resampling schemes differ more.
delta = 0.1 # step size & transition noise (variance)
sigma = 1.0 # observation noise (s.d.)
# generate realisation
mysimobs = ousim(T, delta, sigma, false)

# plot observations (optional: provides reference e.g. for flat portions of tau fn ~ outliers):
scatter(0:T, mysimobs, marker=(:circle, 0.6, :black), leg=false)

# Calculate tau values for all resampling schemes
mysimcNs = smc_example_compare(N, T, mysimobs, ouinit, outransition, oupotential)
sumcN_mn = cumsum(mysimcNs.cN_mn)
sumcN_strat = cumsum(mysimcNs.cN_strat)
sumcN_syst = cumsum(mysimcNs.cN_syst)
sumcN_resmn = cumsum(mysimcNs.cN_resmn)
sumcN_resstrat = cumsum(mysimcNs.cN_resstrat)
sumcN_ressyst = cumsum(mysimcNs.cN_ressyst)

# Plot concurrently tau_N(t) vs t for all resampling schemes
plot([0;sumcN_mn], 0:T, line=(:steppre, 0.6, :black), lab="multinomial", legend=:bottomright, xaxis="t", yaxis="tau_N(t)")
plot!([0;sumcN_strat], 0:T, line=(:steppre, 0.6, :green), lab="stratified")
plot!([0;sumcN_syst], 0:T, line=(:steppre, 0.6, :blue), lab="systematic")
plot!([0;sumcN_resmn], 0:T, line=(:steppre, 0.6, :red), lab="residual-mn")
plot!([0;sumcN_resstrat], 0:T, line=(:steppre, 0.6, :purple), lab="residual-strat")
plot!([0;sumcN_ressyst], 0:T, line=(:steppre, 0.6, :orange), lab="residual-syst")
title!("T = $(T), N = $(N), delta = $(delta), sigma = $(sigma)")

# Export plot as PDF
# savefig("mylovelyplot.pdf")
