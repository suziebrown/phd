w_vals = 0:0.01:1 # values of w1 to evaluate at
nvals = length(w_vals)

myprod = Array{Float64, 2}(undef, nvals, nvals)

for i in 1:nvals
    w1 = w_vals[i]
    for j in 1:nvals
        w2 = w_vals[j]
        w3 = 1 - w1 - w2
        myprod[i,j] = w1*w2*w3
    end
end

using LinearAlgebra
myprod = myprod[nvals:-1:1 , :]
myprod = LowerTriangular(myprod)

using Plots
heatmap(w_vals, w_vals, myprod, xaxis=false)
contour(w_vals, w_vals, myprod, xaxis=false)
