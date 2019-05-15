w = 0:0.001:1
N =20
plot(w, floor.(N.*w), leg=false)
plot!(w, N .* w.^2)

fw = w.^2 .*(N*(2*N-1)) + floor.(N.*w) .* (1 .+ 2 .* floor.(N.*w) .- (4*N) .*w)
plot(w, fw)

fw2 = (N .* w - floor.(N .* w)).^2 - N .* w.^2 + floor.(N .* w)
plot(w, fw2)
plot!(w, (-w.^2 .+ w .- ((N-1)/N^2)).*N)

plot(w, fw2 .- (-w.^2 .+ w .- ((N-1)/N^2)).*N)
plot(w, 2 .* w.^2 -w)
