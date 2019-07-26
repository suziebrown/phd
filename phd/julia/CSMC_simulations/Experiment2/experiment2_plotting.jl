# --- Plot results from Experiment 2 ---

N = 1024
nvals =  [2, 4, 8, 16, 32, 64, 128, 256, 512, 1023]
nrep = 100

mycolours = [:royalblue3, :midnightblue, :indigo, :purple, :maroon4, :maroon]
mymarkers = [:circle, :xcross, :cross, :diamond, :hexagon]

mean4exc = [0.0879199, 0.537529, 0.219902, 1.24404, 1.10611, 0.526982, 0.630977, 0.971172, 2.67925, 1.73682]
l4exc = [0.0078125, 0.0517578, 0.0683105, 0.0692383, 0.0818848, 0.0751953, 0.0800781, 0.0839844, 0.0849609, 0.0828613]
u4exc = [0.178076, 0.289941, 0.319531, 2.58462, 0.319678, 24.7269, 0.146533, 0.157031, 0.231396, 0.256152]
noob4exc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
mean4inc = [21.3993, 21.4392, 23.0428, 29.979, 22.5569, 20.9088, 24.5245, 21.3731, 24.8118, 20.3262]
l4inc = [1.35229, 1.86772, 1.60059, 1.59985, 1.23042, 1.62314, 2.24395, 2.58452, 1.10347, 1.05688]
u4inc =  [5.88965, 34.7422, 64.1732, 32.5115, 35.9698, 25.7218, 10.6246, 9.53301, 14.8697, 118.405]
noob4inc = [0, 0, 0, 1, 2, 0, 0, 0, 0, 0]
ratio4 = mean4inc./mean4exc

mean3exc = [0.192188, 0.451719, 0.535117, 0.492539, 0.894053, 0.799824, 0.732559, 0.935713, 1.33485, 0.786318]
l3exc = [0.00385742, 0.0504883, 0.0631836, 0.0739258, 0.080957, 0.0820313, 0.0859375, 0.0859375, 0.0720703, 0.0820313]
u3exc = [0.147754, 2.37998, 0.122412, 0.158398, 4.58882, 0.429443, 2.79771, 10.2312, 0.694824, 1.78765]
noob3exc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
mean3inc = [1.55482, 1.90914, 2.13989, 1.89823, 2.26937, 1.98381, 2.30793, 1.83642, 2.34462, 1.87897]
l3inc = [0.152832, 0.254736, 0.180029, 0.244971, 0.240869, 0.279248, 0.291553, 0.231396, 0.252686, 0.233594]
u3inc =  [2.5167, 2.64531, 5.91934, 4.76714, 4.82373, 3.98267, 2.47788, 3.74448, 4.77314, 4.48135]
noob3inc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ratio3 = mean3inc./mean3exc

mean2exc = [0.146289, 0.292432, 0.28292, 0.319795, 0.392109, 0.360156, 0.375117, 0.383184, 0.455303, 0.396387]
l2exc = [0.0124512, 0.0497559, 0.0663086, 0.0829102, 0.096582, 0.0897949, 0.109375, 0.0976074, 0.11499, 0.100488]
u2exc = [0.115234, 0.930078, 0.533691, 0.367334, 0.244971, 0.596094, 0.527637, 0.17832, 0.308057, 0.587012]
noob2exc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
mean2inc = [0.279463, 0.342441, 0.403711, 0.388857, 0.413027, 0.392373, 0.378662, 0.393125, 0.455303, 0.400986]
l2inc = [0.0435059, 0.0722656, 0.114551, 0.0976563, 0.104199, 0.0976563, 0.128857, 0.101318, 0.11499, 0.105225]
u2inc =  [0.221143, 0.936816, 2.18657, 0.458936, 0.31084, 0.639355, 0.527637, 0.17832, 0.308057, 0.587012]
noob2inc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ratio2 = mean2inc./mean2exc

mean1exc = [0.0976172, 0.136104, 0.155693, 0.162002, 0.189834, 0.19126, 0.16917, 0.1775, 0.177393, 0.178838]
l1exc = [0.0125488, 0.038916, 0.0623047, 0.0808594, 0.0691895, 0.0820313, 0.0692871, 0.0818848, 0.0839844, 0.0858887]
u1exc = [0.426416, 0.3896, 0.238525, 0.294336, 0.330957, 0.179932, 0.27124, 0.161914, 0.314697, 0.337988]
noob1exc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
mean1inc = [0.0974805, 0.142334, 0.151338, 0.159648, 0.189795, 0.191328, 0.16917, 0.1775, 0.177393, 0.178838]
l1inc = [0.0097168, 0.0398926, 0.0644043, 0.0770996, 0.0691895, 0.0820313, 0.0692871, 0.0818848, 0.0839844, 0.0858887]
u1inc =  [0.070166, 0.174414, 0.238525, 0.298047, 0.330957, 0.179932, 0.27124, 0.161914, 0.314697, 0.337988]
noob1inc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ratio1 = mean1inc./mean1exc

# with ribbons
plot()
plot!(nvals, mean4exc, ribbon=(l4exc, u4exc), label="MAP+4SD:exc", marker=mymarkers[5], markerstrokewidth=0, seriescolor=5,
    xaxis=:log10, title="including/excluding immortal particle: N=$N, nrep=$nrep",
    xlabel="n", ylabel="tree height /N", leg=:topleft)
plot!(nvals[noob4inc.==0], mean4inc[noob4inc.==0], ribbon=(l4inc[noob4inc.==0], u4inc[noob4inc.==0]), label="MAP+4SD:inc", marker=mymarkers[5], markerstrokewidth=0, seriescolor=5, line=:dash)
plot!(nvals, mean3exc, ribbon=(l3exc, u3exc), label="MAP+3SD:exc", marker=mymarkers[4], markerstrokewidth=0, seriescolor=4)
plot!(nvals, mean3inc, ribbon=(l3inc, u3inc), label="MAP+3SD:inc", marker=mymarkers[4], markerstrokewidth=0, seriescolor=4, line=:dash)
plot!(nvals, mean2exc, ribbon=(l2exc, u2exc), label="MAP+2SD:exc", marker=mymarkers[3], markerstrokewidth=0, seriescolor=3)
plot!(nvals, mean2inc, ribbon=(l2inc, u2inc), label="MAP+2SD:inc", marker=mymarkers[3], markerstrokewidth=0, seriescolor=3, line=:dash)
plot!(nvals, mean1exc, ribbon=(l1exc, u1exc), label="MAP+1SD:exc", marker=mymarkers[2], markerstrokewidth=0, seriescolor=2)
plot!(nvals, mean1inc, ribbon=(l1inc, u1inc), label="MAP+1SD:inc", marker=mymarkers[2], markerstrokewidth=0, seriescolor=2, line=:dash)

# separately for each nsd
plot()
plot!(nvals, mean4exc, ribbon=(l4exc, u4exc), label="MAP+4SD:exc", marker=mymarkers[5], markerstrokewidth=0, seriescolor=5,
    xaxis=:log10, title="including/excluding immortal particle: N=$N, nrep=$nrep, nsd=4", titlefontsize=10,
    xlabel="n", ylabel="tree height /N", leg=:topleft)
plot!(nvals[noob4inc.==0], mean4inc[noob4inc.==0], ribbon=(l4inc[noob4inc.==0], u4inc[noob4inc.==0]), label="MAP+4SD:inc", marker=mymarkers[5], markerstrokewidth=0, seriescolor=5, line=:dash)

plot()
plot!(nvals, mean3exc, ribbon=(l3exc, u3exc), label="MAP+3SD:exc", marker=mymarkers[4], markerstrokewidth=0, seriescolor=4,
    xaxis=:log10, title="including/excluding immortal particle: N=$N, nrep=$nrep, nsd=3", titlefontsize=10,
    xlabel="n", ylabel="tree height /N", leg=:topleft)
plot!(nvals, mean3inc, ribbon=(l3inc, u3inc), label="MAP+3SD:inc", marker=mymarkers[4], markerstrokewidth=0, seriescolor=4, line=:dash)

plot()
plot!(nvals, mean2exc, ribbon=(l2exc, u2exc), label="MAP+2SD:exc", marker=mymarkers[3], markerstrokewidth=0, seriescolor=3,
    xaxis=:log10, title="including/excluding immortal particle: N=$N, nrep=$nrep, nsd=2", titlefontsize=10,
    xlabel="n", ylabel="tree height /N", leg=:topleft)
plot!(nvals, mean2inc, ribbon=(l2inc, u2inc), label="MAP+2SD:inc", marker=mymarkers[3], markerstrokewidth=0, seriescolor=3, line=:dash)

plot()
plot!(nvals, mean1exc, ribbon=(l1exc, u1exc), label="MAP+1SD:exc", marker=mymarkers[2], markerstrokewidth=0, seriescolor=2,
    xaxis=:log10, title="including/excluding immortal particle: N=$N, nrep=$nrep, nsd=1", titlefontsize=10,
    xlabel="n", ylabel="tree height /N", leg=:topleft)
plot!(nvals, mean1inc, ribbon=(l1inc, u1inc), label="MAP+1SD:inc", marker=mymarkers[2], markerstrokewidth=0, seriescolor=2, line=:dash)


# without ribbons
plot()
plot!(nvals, mean4exc, label="MAP+4SD:exc", marker=mymarkers[5], markerstrokewidth=0, seriescolor=5,
    xaxis=:log10, yaxis=:log10, title="including/excluding immortal particle: N=$N, nrep=$nrep", titlefontsize=10,
    xlabel="n", ylabel="tree height /N", xlims=(1,10^4.5), leg=:topright)
plot!(nvals[noob4inc.==0], mean4inc[noob4inc.==0], label="MAP+4SD:inc", marker=mymarkers[5], markerstrokewidth=0, seriescolor=5, line=:dash)
plot!(nvals, mean3exc, label="MAP+3SD:exc", marker=mymarkers[4], markerstrokewidth=0, seriescolor=4)
plot!(nvals, mean3inc, label="MAP+3SD:inc", marker=mymarkers[4], markerstrokewidth=0, seriescolor=4, line=:dash)
plot!(nvals, mean2exc, label="MAP+2SD:exc", marker=mymarkers[3], markerstrokewidth=0, seriescolor=3)
plot!(nvals, mean2inc, label="MAP+2SD:inc", marker=mymarkers[3], markerstrokewidth=0, seriescolor=3, line=:dash)
plot!(nvals, mean1exc, label="MAP+1SD:exc", marker=mymarkers[2], markerstrokewidth=0, seriescolor=2)
plot!(nvals, mean1inc, label="MAP+1SD:inc", marker=mymarkers[2], markerstrokewidth=0, seriescolor=2, line=:dash)

# plot ration of include:exclude
plot()
plot!(nvals[noob4inc.==0], ratio4[noob4inc.==0], label="MAP+4SD", marker=mymarkers[5], markerstrokewidth=0, seriescolor=5,
    xaxis=:log10, yaxis=:log10, title="ratio including:excluding immortal particle: N=$N, nrep=$nrep", titlefontsize=10,
    xlabel="n", ylabel="ratio of tree heights", leg=:topright)
plot!(nvals, ratio3, label="MAP+3SD", marker=mymarkers[4], markerstrokewidth=0, seriescolor=4)
plot!(nvals, ratio2, label="MAP+2SD", marker=mymarkers[3], markerstrokewidth=0, seriescolor=3)
plot!(nvals, ratio1, label="MAP+1SD", marker=mymarkers[2], markerstrokewidth=0, seriescolor=2)

#savefig("CSMC_treeheight_incexc_nsd4.pdf")
