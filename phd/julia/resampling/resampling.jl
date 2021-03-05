### Implementation of different resampling schemes
# OBSOLETE: a complete resampling implementation can be found at:
# phd/julia/resampling/resampling_CP2020

using Distributions

N=4

# Randomly generate the vector of weights
function getw(N)
    w = rand(Float64, N)
    sumw = sum(w)
    w/sumw
end


# Multinomial resampling
function resam_mn(N,w)
    rand(Multinomial(N,w))
end


# Residual-multinomial resampling
function resam_resmn(N,w)
    floors = Int.(floor.(N.*w))
    resids = (N.*w - floors)
    R = sum(resids)
    resids = resids/R
    floors + rand(Multinomial(Int(R),resids))
end


# Stratified resampling


# Systematic resampling
