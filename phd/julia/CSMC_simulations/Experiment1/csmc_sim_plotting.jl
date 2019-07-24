#---- Plot results of csmc simulations ----

# zerosd_mean = [625.24, 1030.3, 1109.64, 1231.71, 1286.67, 1256.7, 1353.16, 1289.33, 1326.65, 1389.62, 1274.49, 1378.2, 1293.91]
# zerosd_lower = [30.1, 252.8, 333.7, 443.7, 602.4, 538.75, 639.6, 587.4, 589.95, 559.25, 635.35, 556.9, 579.4]
# zerosd_upper = [1278.35, 2555.4, 1321.25, 1638.1, 1379.75, 1823.9, 1620.55, 833.7, 2826.45, 1883.6, 946.45, 2656.2, 2315.15]
#
# onesd_mean = [890.74, 1274.41, 1623.46, 1738.08, 1754.06, 1800.82, 1536.5, 1672.95, 1649.34, 1989.71, 1762.78, 1828.11, 1787.68]
# onesd_lower = [80.75, 245.85, 434.4, 624.25, 606.75, 693.8, 684.25, 689.3, 642.2, 712.7, 641.55, 670.75, 692.75]
# onesd_upper = [1058.55, 2281.05, 2831.15, 2378.95, 4505.75, 1903.45, 2497.0, 1578.15, 3216.95, 2753.1, 1603.15, 2689.7, 1664.4]
#
# twosd_mean = [1704.7, 2235.71, 3098.44, 2804.35, 3224.21, 3508.56, 3425.71, 3550.22, 3780.95, 3718.0, 3309.59, 3377.76, 3490.16]
# twosd_lower = [49.75, 248.85, 555.7, 691.3, 779.45, 751.7, 917.5, 1121.95, 1288.8, 1100.95, 952.7, 879.5, 993.6]
# twosd_upper = [6501.7, 4870.7, 6067.5, 5740.45, 2498.9, 5242.55, 2711.8, 6773.9, 4272.65, 5121.35, 8000.95, 3916.85, 4138.85]
#
# threesd_mean = [1303.58, 2021.69, 2930.37, 3351.76, 3235.06, 3929.04, 4223.6, 4162.45, 4256.62, 5497.46, 6372.26, 7230.03, 7683.48]
# threesd_lower = [18.0, 259.95, 347.8, 608.0, 746.5, 670.7, 610.35, 669.75, 789.85, 930.0, 855.15, 985.9, 1214.35]
# threesd_upper = [1355.55, 2694.5, 9655.4, 9915.95, 9743.55, 10000.0, 5649.6, 10000.0, 9587.5, 10000.0, 9088.95, 10000.0, 10000.0]
#
# foursd_mean = [776.73, 1319.05, 1331.11, 1497.77, 1692.49, 1587.48, 1784.62, 2281.51, 2514.92, 2517.93, 3880.15, 6241.81, 7287.09]
# foursd_lower = [40.85, 366.3, 487.95, 636.85, 620.25, 540.6, 673.65, 723.0, 690.0, 652.5, 706.6, 779.4, 816.65]
# foursd_upper = [1188.8, 2924.4, 1397.35, 3980.2, 3036.6, 2503.3, 969.1, 10000.0, 9604.25, 9648.05, 10000.0, 10000.0, 10000.0]
#
# plot(nvals, zerosd_mean/N, ribbon=((zerosd_lower)/N,(zerosd_upper)/N), label="MAP", fill=:red, fillalpha=0.25, leg=:topleft, xaxis=:log10, line=(:red), marker=(:red), markerstrokecolor=:red, title="CSMC treeheight, N=8192, nrep=100", xlabel="n", ylabel="tree height /N")
# plot!(nvals, onesd_mean/N, ribbon=((onesd_lower)/N,(onesd_upper)/N), label="MAP+1SD", fill=:purple, fillalpha=0.25, line=(:purple), marker=(:purple), markerstrokecolor=:purple)
# plot!(nvals, twosd_mean/N, ribbon=((twosd_lower)/N,(twosd_upper)/N), label="MAP+2SD", fill=:blue, fillalpha=0.25, line=(:blue), marker=(:blue), markerstrokecolor=:blue)
# plot!(nvals, threesd_mean/N, ribbon=((threesd_lower)/N,(threesd_upper)/N), label="MAP+3SD", fill=:green, fillalpha=0.25, line=(:green), marker=(:green), markerstrokecolor=:green)
# plot!(nvals, foursd_mean/N, ribbon=((foursd_lower)/N,(foursd_upper)/N), label="MAP+4SD", fill=:yellow, fillalpha=0.25, line=(:yellow), marker=(:yellow), markerstrokecolor=:yellow)

#savefig("CSMC_treeheight_A.pdf")
#
# N=8192
#
# mycolours = [:royalblue3, :midnightblue, :indigo, :purple, :maroon4, :maroon]
#
# zerosd = [726.102, 1016.82, 1155.61, 1201.44, 1362.67, 1348.3, 1323.37, 1340.24, 1359.86, 1342.16, 1338.39, 1307.77, 1321.49]
# onesd = [885.014, 1288.91, 1581.94, 1679.28, 1718.06, 1643.04, 1649.05, 1757.22, 1762.51, 1759.39, 1671.97, 1703.9, 1750.46]
# twosd = [1451.78, 2279.75, 2785.59, 3167.54, 3721.42, 3295.2, 3471.46, 3606.28, 3763.73, 3657.91, 3731.39, 3689.7, 3799.12]
# threesd = [1331.69, 1959.96, 2967.38, 3164.86, 3278.85, 4246.1, 4097.17, 4489.82, 4871.32, 5202.34, 5850.78, 6302.09, 7287.7]
# foursd = [909.976, 1357.13, 1751.93, 1808.43, 1791.75, 1964.09, 2068.96, 2259.14, 2649.8, 2820.46, 3940.88, 5179.25, 7102.77]
# fivesd = [814.142, 1205.01, 1480.15, 1516.35, 1533.83, 1742.18, 1667.74, 1929.58, 2259.44, 2671.45, 3427.94, 5107.16, 6883.01]
#
# plot()
# plot!(nvals, zerosd/N, label="MAP", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[1], title="CSMC treeheight, N=8192", xlabel="n", ylabel="mean tree height /N")
# plot!(nvals, onesd/N, label="MAP+1SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[2])
# plot!(nvals, twosd/N, label="MAP+2SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[3])
# plot!(nvals, threesd/N, label="MAP+3SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[4])
# plot!(nvals, foursd/N, label="MAP+4SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[5])
# plot!(nvals, fivesd/N, label="MAP+5SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[6])

#savefig("CSMC_treeheight_n500_B.pdf")


N=1024
nvals =  [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

mycolours = [:royalblue3, :midnightblue, :indigo, :purple, :maroon4, :maroon]

zerosd = [0.093502, 0.136135, 0.158344, 0.17349, 0.180461, 0.178074, 0.184531, 0.181078, 0.183539, 0.19099]
zerosd_upper = [0.194287, 0.182324, 0.209668, 0.211279, 0.256592, 0.449951, 0.301025, 0.153857, 0.371729, 0.287158]
zerosd_lower = [0.0078125, 0.0361328, 0.0634277, 0.0673828, 0.0722656, 0.0751465, 0.0780762, 0.0712891, 0.0712891, 0.0799805]
#onesd = [885.014, 1288.91, 1581.94, 1679.28, 1718.06, 1643.04, 1649.05, 1757.22, 1762.51, 1759.39, 1671.97, 1703.9, 1750.46]
twosd = [0.19666, 0.328717, 0.43217, 0.42248, 0.481, 0.48666, 0.508152, 0.530256, 0.526613, 0.505318]
twosd_upper = [0.150684, 0.105713, 0.280176, 0.432275, 1.4625, 0.626221, 0.840869, 0.518994, 0.526074, 0.758936]
twosd_lower = [0.00976563, 0.0478516, 0.0712402, 0.0956055, 0.125928, 0.135596, 0.157178, 0.15625, 0.159033, 0.159814]
#
threesd = [0.322908, 0.429732, 0.620645, 0.637828, 0.909016, 0.91867, 1.07224, 1.47865, 1.83172, 2.3112]
threesd_upper = [0.0724609, 0.828809, 0.277588, 0.257178, 3.16104, 4.49238, 0.25083, 2.75547, 0.468555, 2.3708]
threesd_lower = [0.00878906, 0.046875, 0.0683594, 0.0761719, 0.0868164, 0.0908203, 0.0996094, 0.120898, 0.151172, 0.322021]
#
foursd = [0.295217, 0.482889, 1.03172, 1.58346, 1.44648, 2.83352, 4.35965, 6.76921, 12.8193, 24.7134]
foursd_upper = [0.157666, 0.323926, 0.210889, 0.437744, 0.331445, 0.189062, 0.227393, 8.1377, 30.2178, 26.8226]
foursd_lower = [0.00976563, 0.0428711, 0.0673828, 0.0742188, 0.0849121, 0.0986328, 0.0839355, 0.0888672, 0.105273, 1.94282]
#fivesd = [814.142, 1205.01, 1480.15, 1516.35, 1533.83, 1742.18, 1667.74, 1929.58, 2259.44, 2671.45, 3427.94, 5107.16, 6883.01]

#extra lines
plot()
plot!(nvals, zerosd, label="MAP", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[1], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")
plot!(nvals, onesd, label="MAP+1SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[2])
plot!(nvals, twosd, label="MAP+2SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[3])
plot!(nvals, threesd, label="MAP+3SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[4])
plot!(nvals, foursd, label="MAP+4SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[5])
plot!(nvals, fivesd, label="MAP+5SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[6])

#without ribbons
plot()
plot!(nvals, zerosd, label="MAP", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[1], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")
plot!(nvals[1:9], twosd[1:9], label="MAP+2SD", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[3], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")
plot!(nvals, threesd, label="MAP+3SD", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[4], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")
plot!(nvals[1:9], foursd[1:9], label="MAP+4SD", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[5], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")

# with ribbons
plot()
plot!(nvals, zerosd, ribbon=(zerosd_lower,zerosd_upper), label="MAP", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[1], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")
plot!(nvals[1:9], twosd[1:9], ribbon=(twosd_lower[1:9],twosd_upper[1:9]), label="MAP+2SD", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[3], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")
plot!(nvals[1:9], foursd[1:9], ribbon=(foursd_lower[1:9],foursd_upper[1:9]), label="MAP+4SD", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[5], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")
plot!(nvals, threesd, ribbon=(threesd_lower,threesd_upper), label="MAP+3SD", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[4], title="CSMC treeheight, N=1024", xlabel="n", ylabel="mean tree height /N")

#savefig("CSMC_treeheight_n500_B.pdf")
