# preliminary simulations (using full matrix storage of ancestry)
#
using StatsBase, Random, Distributions, Plots

##--- tree generation

# set constants
delta = 0.04 # noise variance in AR(1) process
sigma = 0.2 # noise s.d. in observations
T = Int64(100) # number of generations/time steps
N = UInt16(20) # total number of particles

# generate observations & immortal trajectory
observations = ousim(T, delta, sigma, false)
immpos = ourts(delta, sigma, observations).mean

# run csmc
csmcout = csmc_fullstore(N, T, observations, ouinit, outransition, oupotential, immpos)
immleaf = csmcout.immortal[T+1]

##--- subtree sampling

# set constants
n = 2 # number of leaves in sampled subtree
nrep = 10

height = Array{Int64, 1}(undef, nrep)
incimm = Array{Bool, 1}(undef, nrep)

for i in 1:nrep
    leaves = sample(1:N, n)
    incimm[i] = immleaf in leaves
    height[i] = mrca_fullstore(csmcout.parents, leaves)
end

println(height)
println(incimm)
meanimm = mean(height[incimm])
meannot = mean(height[.!incimm])
