### ILLUSTRATIVE PLOTS SHOWING ACENSTRAL DEGENERACY
#   also comparing this between multi, adaptive multi, lowvar, adaptive lowvar.
#   ALSO LOAD ALL FUNCTIONS FROM:
#       phd/julia/OU_model.jl
#       phd/julia/resampling/resampling_CP2020.jl
using Random, Distributions, StatsBase, Plots, Statistics


function smc_example_anc(N::Int64, T::Int64, observations::Array{Float64,1}, initialsam::Function, transition::Function, potential::Function)
    ## pre-allocate memory
    parents = Array{Int64, 2}(undef, T, N)
    positions = Array{Float64, 1}(undef, N)
    weights = Array{Float64, 1}(undef, N)

    ## initialise
    positions = initialsam(N, delta, sigma)
    weights = potential(positions, observations[1], delta, sigma)
    weights = weights / sum(weights) ## normalise weights

    ## pre-calculate random numbers for OU transitions
    rands = reshape(rand(Normal(), N*T), T, N)

    for t in 1:T # note: index t+1 corresponds to generation t
        ## choose parents (uncomment if statement for adpative resampling)
        #ess = (sum(weights.^2))^(-1)
        #if ess > N/2
            parents[t, :] = resam_mn(N, weights) ## change resampling scheme as desired
        #else
        #    parents[t, :] = collect(1:N)
        #end
        ## update positions
        positions = outransition(positions[parents[t, :]], delta, sigma, rands[t, :])
        ## compute weights
        weights = potential(positions, observations[t+1], delta, sigma)
        weights = weights / sum(weights)
    end

    return parents
end

## now with all rands set outside of function, for comparisons:
function smc_anc_compare(N::Int64, T::Int64, observations::Array{Float64,1}, transition::Function, potential::Function, initpos::Array{Float64,1}, rands::Array{Float64,2})
    ## pre-allocate memory
    parents = Array{Int64, 2}(undef, T, N)
    weights = Array{Float64, 1}(undef, N)
    ## initialise
    positions = initpos
    weights = potential(positions, observations[1], delta, sigma)
    weights = weights / sum(weights) ## normalise weights

    for t in 1:T # note: index t+1 corresponds to generation t
        ## choose parents (uncomment if statement for adpative resampling)
        #ess = (sum(weights.^2))^(-1)
        #if ess > N/2
            parents[t, :] = resam_syst(N, weights) ## change resampling scheme as desired
        #else
        #    parents[t, :] = collect(1:N)
        #end
        ## update positions
        positions = outransition(positions[parents[t, :]], delta, sigma, rands[t, :])
        ## compute weights
        weights = potential(positions, observations[t+1], delta, sigma)
        weights = weights / sum(weights)
    end

    return parents
end

function trace_ancestry(N::Int64, T::Int64, parents::Array{Int64,2})
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



### The comparison

# set-up
T=50 #time horizon
N=25 # number of particles
delta = 0.1 # step size & transition noise (variance)
sigma = 0.5 # observation noise (s.d.)
mysimobs = ousim(T, delta, sigma, false)
rands = reshape(rand(Normal(), N*T), T, N)
initpos = ouinit(N, delta, sigma)

# run the thing (edit fn & name parents var differently for each thing to compare)
myparents = smc_anc_compare(N, T, mysimobs, outransition, oupotential, initpos, rands)
myanc = trace_ancestry(N,T,myparents)

# plot results
myancplot = plot(1:T+1, myanc[:,1], line=(:black), legend=false, xaxis="t", yaxis="ancestral index", ylims=(-0.5,N+0.5))
for i in 2:N
    plot!(1:T+1, myanc[:,i], line=:black)
end
title!("ancestral degeneracy") ## change plot title for different cases
myancplot
#savefig("mylovelyplot.pdf") ## change file name for different cases
