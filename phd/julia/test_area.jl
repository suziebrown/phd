# TEST AREA
#
using StatsBase, Random, Distributions

delta = 1.0
sigma = 0.1
T = 10
N = 5

testobs = ousim(T, delta, sigma, false)
immpos = ousim(T, delta, sigma, false)

testres = csmc_fullstore(N, T, testobs, ouinit, outransition, oupotential, immpos)

# sub-tree height (MRCA) searching function
mrca_fullstore(testres.parents, [1,2,4])
