# Plots from Experiment 3

nvals = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1023]

# proportion of samples finding the immortal line
prop0 = [78.0, 92.0, 97.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0] /100
prop1 = [67.0, 92.0, 96.0, 99.0, 99.0, 100.0, 100.0, 100.0, 100.0, 100.0] /100
prop2 = [38.0, 77.0, 72.0, 85.0, 88.0, 94.0, 95.0, 96.0, 99.0, 99.0] /100
prop3 = [5.0, 16.0, 29.0, 28.0, 38.0, 39.0, 45.0, 50.0, 52.0, 48.0] /100
prop4 = [0.0, 1.0, 0.0, 5.0, 2.0, 2.0, 4.0, 3.0, 2.0, 3.0] /100

plot()
plot!(nvals, prop0, label="MAP", xaxis=:log10, xlabel="n", ylabel="proportion of samples finding immortal line")
plot!(nvals, prop1, label="MAP+1SD")
plot!(nvals, prop2, label="MAP+2SD")
plot!(nvals, prop3, label="MAP+3SD")
plot!(nvals, prop4, label="MAP+4SD")
title!("N=1024, nrep=100")

savefig("proportion_coalescing_immortal.pdf")
