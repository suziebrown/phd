#---- c_N ----
w = 0:0.001:1
N =10
plot(w, floor.(N.*w), leg=false)
plot!(w, N .* w.^2)

# wrong function, I think
fw = w.^2 .*(N*(2*N-1)) + floor.(N.*w) .* (1 .+ 2 .* floor.(N.*w) .- (4*N) .*w)
plot(w, fw)

# plot Delta_i
fw2 = (N .* w - floor.(N .* w)).^2 - N .* w.^2 + floor.(N .* w)
plot(w, fw2, leg=false, title="N=10", xlab="w")
plot!(w, (-w.^2 .+ w .- ((N-1)/N^2)).*N)

# distance from fitted parabola
plot(w, fw2 - (-w.^2 .+ w .- ((N-1)/N^2)).*N)
plot(w, 2 .* w.^2 -w)


#---- D_N ----
w = 0:0.001:0.2
N =100
k = floor.(N*w)

# part without crossterms - RR alone and difference Mn-RR
fw3 = (-2) * floor.(N*w).^2 + 3*N* w .* floor.(N*w).^2 - 2* floor.(N*w).^3 + N*w .* floor.(N*w)
plot(w, fw3, leg=false)
plot!(w, N*(N-1)*(N-2)* w.^3 + 2*N*(N-1)* w.^2)

plot(w, N*(N-1)*(N-2)* w.^3 + 2*N*(N-1)* w.^2 - fw3, leg=false)

#part with crossterms - RR alone
fw4 = floor.(N*w) .* (N*(2*N*w .-1) + floor.(N*w) * (-1)*(N+1)+ (floor.(N*w).^2) .* (2 .- 4*N*w) + 3 * (floor.(N*w).^3) )
plot(w, fw4)

# whole thing - UB for RR
fw5a = k*N + k.^2 *(N^2 +1)/N + k.^3 *((3*N-2)/N) - k.^4 * (N-1)/N
plot(w, fw5a, label="RR upper", leg=:bottomleft)
title!("Issue with D_N ordering for some w")
xlabel!("w")
# whole thing - LB for Mn
fw5b = k.^3 * (N-1)*(N-2)/N^2 + k.^2 * 2*(N-1)/N
plot!(w, fw5b, label="Mn lower")
# whole thing - bound on difference
fw5 = k.^4 *(N-1)/N + k.^3 * ((2/N^2) - (1/N) -2) - k.^2 * (N-3)*(N+1)/N - k *N
plot(w, fw5)

#---- save plot ----
#savefig("delta_DN_N10.pdf")
