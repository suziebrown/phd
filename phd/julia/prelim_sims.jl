# preliminary simulations (using full matrix storage of ancestry)
#
using StatsBase, Random, Distributions, Plots

##--- tree generation

# set constants
delta = 0.1 # noise variance in AR(1) process
sigma = 0.1 # noise s.d. in observations
T = Int64(10000) # number of generations/time steps
N = Int64(8192) # total number of particles

# generate observations & immortal trajectory
observations = ousim(T, delta, sigma, false)
immpos = ourts(delta, sigma, observations).mean

# run csmc
#csmcout = csmc_fullstore(N, T, observations, ouinit, outransition, oupotential, immpos)
#immleaf = csmcout.immortal[T+1]

##--- subtree sampling

# set constants
nvals = [2,4,8,16,32,64,128,256,512,1024,2048,4096,8192] # number of leaves in sampled subtree
nrep = 10

# initialise local variables
height = Array{Int64, 1}(undef, nrep)
#incimm = Array{Bool, 1}(undef, nrep)

# initialise output variables
meanall = Array{Float64, 1}(undef, length(nvals))
#meanimm = Array{Float64, 1}(undef, length(nvals))
#meannot = Array{Float64, 1}(undef, length(nvals))
sdall = Array{Float64, 1}(undef, length(nvals))
noob = zeros(Int64, length(nvals))

for j in 1:length(nvals)
    # uniformly sample 'nrep'x subtrees of size 'nvals[j]' & calculate height
    for i in 1:nrep
        # run csmc
        csmcout = csmc_fullstore(N, T, observations, ouinit, outransition, oupotential, immpos)
        # sample a subtree
        leaves = sample(1:N, nvals[j])
        #incimm[i] = immleaf in leaves
        height[i] = mrca_fullstore(csmcout.parents, leaves)

        # error catching
        noob[j] += sum(height==T)
    end
    # mean heights for all samples, samples inc. immortal particle, & samples w/o immortal particle.
    meanall[j] = mean(height)
    #meanimm[j] = mean(height[incimm])
    #meannot[j] = mean(height[.!incimm])
    sdall[j] = std(height)
end

# catching error of T being too small for N, i.e. reports no. of cases where treeheight was T
println("number of cases hitting limit was ", sum(noob))
println("means were: ", meanall)
println("SDs were: ", sdall)
# plot output (ribbon shows +/- 1 standard error)
plot(nvals, meanall/N, ribbon=(sdall*nrep^(-0.5)/N, sdall*nrep^(-0.5)/N), fill=:purple, fillalpha=0.25, leg=false, xaxis=:log10, line=(:purple), marker=(:purple), markerstrokecolor=:purple, title="CSMC treeheight, immortal=MAP, N=$N", xlabel="n", ylabel="average tree height /N")
#plot!(nvals, meanimm, label="inc. immortal")
#plot!(nvals, meannot, label="w/o immortal")



#---- ---- ---- ---- ---- ---- ---- ----
# STORING RESULTS FOR T=10^4, N=2048:
#meanallsafe = meanall
#sdallsafe = sdall
