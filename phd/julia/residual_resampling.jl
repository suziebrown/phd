using Random, Distributions

n = 10 # number of particles
reps = 4 # number of weight vectors to sample
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

v_res_det = Array{Int64, 1}(undef, n)
v_res_rand = Array{Int64, 1}(undef, n)
resid = Array{Float64, 1}(undef, n)

for i in 1:reps
    v_mn[i,:] = rand(Multinomial(n, w[i,:])) # multinomial resampling

    v_res_det = floor.(Int, n * w[i,:]) # deterministic part
    resid = w[i,:] .- (v_res_det /n) # residual weights
    # for j in 1:n
    #     v_res_det[j] = floor(Int, n * w[i,j]) # deterministic part
    #     resid[j] = w[i,j] - (v_res_det[j] /n) # residual weights
    # end
    R = sum(resid)
    resid = resid / R
    R = Int(round(n * R))
    v_res_rand = rand(Multinomial(R, resid)) # random part

    v_res[i,:] = v_res_det + v_res_rand # residual resampling
end
