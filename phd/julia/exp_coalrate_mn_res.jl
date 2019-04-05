### comparing E(c_N^r) and E(c_N^m) for N=2 over varying w_1

w1_vals = 0:0.01:1 # values of w1 to evaluate at
nvals = length(w1_vals)

EcN_mn = Array{Float64, 1}(undef, nvals)
EcN_res = Array{Float64, 1}(undef, nvals)

for i in 1:length(w1_vals)
    # set weights
    w1 = w1_vals[i]
    w2 = 1 - w1
    # probability of having 0,1,2 offspring assigned to parent 1 (multinomial):
    p0_mn = w2 ^2
    p1_mn = 2 * w1 * w2
    p2_mn = w1 ^2
    # probability of having 0,1,2 offspring assigned to parent 1 (residual):
    p0_res = (w1 < 0.5 ? 1 : 0) * ((w2 - 0.5) * 2)
    p1_res = (w1 < 0.5 ? 1 : 0) * (w1 * 2) + (w2 < 0.5 ? 1 : 0) * (w2 * 2)
    p2_res = (w1 > 0.5 ? 1 : 0) * ((w1 - 0.5) * 2)
    # expected coalescence rates
    EcN_mn[i] = p0_mn + p2_mn
    EcN_res[i] = p0_res + p2_res
end

# plot results
using Plots

plot(w1_vals, EcN_mn, lab="multinomial", line=2, legend=:bottomright)
plot!(w1_vals, EcN_res, lab="residual", line=2)
ylabel!("expected coalescence rate")
xlabel!("w1")
title!("dependence of E[c_N] on weights (N=2)")

#savefig("EcN_mn_res_N2.pdf")
