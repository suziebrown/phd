# TEST AREA


##--- test ousim and plotting positions

using Plots
gr()

x = ousim(10, 0.1, 1.0, true)
plot(0:10, x.states, label="states")
plot!(0:10, x.observations, label="observations")

x = ousim(10, 0.1, 1.0, false)
plot(0:10, x, label="observations")


##--- testing csmc algorithm & mrca search algorithm

using StatsBase, Random, Distributions

delta = 1.0
sigma = 0.1
T = 10
N = 5

testobs = ousim(T, delta, sigma, false)
immpos = ousim(T, delta, sigma, false)

testres = csmc_fullstore(N, T, testobs, ouinit, outransition, oupotential, immpos)

mrca_fullstore(testres.parents, [1,2,4])


##--- testing kalman & rts functions

using Plots

delta = 0.1 # noise in AR(1) process
sigma = 1.0 # noise in observations
T = 10

testobs = ousim(T, delta, sigma, true)
kalman = oukalman(delta, sigma, testobs.observations)
rts = ourts(delta, sigma, testobs.observations)

plot(0:T, testobs.states, label="states")
plot!(0:T, testobs.observations, label="obs")
plot!(0:T, kalman.mean, label="kalman")
plot!(0:T, rts.mean, label="rts")
