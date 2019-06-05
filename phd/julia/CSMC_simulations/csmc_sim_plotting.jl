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

N=8192

mycolours = [:royalblue3, :midnightblue, :indigo, :purple, :maroon4]

zerosd = [726.102, 1016.82, 1155.61, 1201.44, 1362.67, 1348.3, 1323.37, 1340.24, 1359.86, 1342.16, 1338.39, 1307.77, 1321.49]
onesd = [885.014, 1288.91, 1581.94, 1679.28, 1718.06, 1643.04, 1649.05, 1757.22, 1762.51, 1759.39, 1671.97, 1703.9, 1750.46]
twosd = [1451.78, 2279.75, 2785.59, 3167.54, 3721.42, 3295.2, 3471.46, 3606.28, 3763.73, 3657.91, 3731.39, 3689.7, 3799.12]
threesd = [1331.69, 1959.96, 2967.38, 3164.86, 3278.85, 4246.1, 4097.17, 4489.82, 4871.32, 5202.34, 5850.78, 6302.09, 7287.7]
foursd = [909.976, 1357.13, 1751.93, 1808.43, 1791.75, 1964.09, 2068.96, 2259.14, 2649.8, 2820.46, 3940.88, 5179.25, 7102.77]

plot()
plot!(nvals, zerosd_mean/N, label="MAP", marker=:auto, markerstrokewidth=0, leg=:topleft, xaxis=:log10, seriescolor=mycolours[1], title="CSMC treeheight, N=8192", xlabel="n", ylabel="mean tree height /N")
plot!(nvals, onesd_mean/N, label="MAP+1SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[2])
plot!(nvals, twosd_mean/N, label="MAP+2SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[3])
plot!(nvals, threesd_mean/N, label="MAP+3SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[4])
plot!(nvals, foursd_mean/N, label="MAP+4SD", marker=:auto, markerstrokewidth=0, seriescolor=mycolours[5])

savefig("CSMC_treeheight_n500.pdf")
