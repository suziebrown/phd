using Random, Distributions

n = 10 # number of particles
reps = 100 # number of weight vectors to sample
w = Array{Float64, 2}(undef, reps, n)

# generate weight vectors
for i in 1:reps
    w[i, :] = rand(Float64, n)
end
# normalise weight vectors
W = sum(w, dims=2).^(-1)
for i in 1:reps
    w[i, :] = w[i,:] * W[i]
end

v_mn = Array{Int64, 2}(undef, reps, n)
v_res = Array{Int64, 2}(undef, reps, n)
n_random = Array{Int64, 1}(undef, reps)
mad_mn = Array{Float64, 1}(undef, reps)
mad_res = Array{Float64, 1}(undef, reps)

# simulate resampling with multinomial & residual
for i in 1:reps
    v_mn[i,:] = rand(Multinomial(n, w[i,:])) # multinomial resampling
    mad_mn[i] = mean(abs.(v_mn[i,:] .- n * w[i,:])) # mean absolute deviation from mean offspring numbers

    # residual resmapling procedure...
    v_res_det = floor.(Int, n * w[i,:]) # deterministic part
    resid = w[i,:] .- (v_res_det /n) # residual weights
    R = sum(resid)
    resid = resid / R
    R = Int(round(n * R))
    n_random[i] = R # how many offspring were assigned randomly
    v_res_rand = rand(Multinomial(R, resid)) # random part

    v_res[i,:] = v_res_det + v_res_rand # residual resampling
    mad_res[i] = mean(abs.(v_res[i,:] .- n * w[i,:]))
end

# compare offspring variance etc...
pairs_mn = Array{Int64, 1}(undef, reps)
binom_mn = Array{Int64,1}(undef, n)
pairs_res = Array{Int64, 1}(undef, reps)
binom_res = Array{Int64,1}(undef, n)

for i in 1:reps
    for j in 1:n
        binom_mn[j] = binomial(v_mn[i,j],2)
        binom_res[j] = binomial(v_res[i,j],2)
    end
    pairs_mn[i] = sum(binom_mn)
    pairs_res[i] = sum(binom_res)
end

n_greater = sum(pairs_res .> pairs_mn)
n_equal = sum(pairs_res .== pairs_mn)
avg_random = mean(n_random)
mad_greater = sum(mad_res .> mad_mn)

println()
println("On average, residual assigned ", avg_random, " out of ", n, " offspring randomly.")
println("Residual merged more pairs than multinomial in ", n_greater, " out of ", reps, " cases.")
println("Residual merged the same number of pairs as multinomial in ", n_equal, " out of ", reps, " cases.")
println("Residual gave higher MAD from expected offspring numbers in ", mad_greater, " out of ", reps, " cases.")
println()
