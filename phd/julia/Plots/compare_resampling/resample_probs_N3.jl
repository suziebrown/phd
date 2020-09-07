### Ternary plots comparing resampling schemes on fixed weights

function prob111mn(x,y)
    6 * x .* y .* (1 .-x .-y)
end

# function prob111res(x::Float64, y::Float64) # scalar inputs only in this version
#     N = 3
#     w = [x; y; (1 -x -y)]
#     flnw = floor.(N * w)
#     nubar = 1 .- flnw
#     wbar = N * w .- flnw
#     R::Int64 = sum(wbar)
#     if (R>0)
#         wbar = wbar/R
#     end
#     return(sum(nubar .< 0) == 0 ? factorial(R) * prod(wbar .^ nubar) : 0)
# end

# attempt at vectorised version (not currently working)
function prob111res(x, y) # set input types
    N = 3
    out = Array{Float64, 1}(undef, length(x))
    w = [x y (1 .-x .-y)]
    flnw = floor.(N * w)
    nubar = 1 .- flnw
    wbar = N * w .- flnw
    R = convert.(Int64, round.(sum(wbar,dims=2)))
    for i in 1:length(x)
        wbar[i,:] = (R[i]>0) ? (wbar[i,:] / R[i]) : (wbar[i,:] * 0)
        out[i] = (sum(nubar[i,:] .< 0) == 0) ? factorial(R[i]) * prod(wbar[i,:] .^ nubar[i,:]) : 0
    end
    return(out)
    #return((sum(nubar .< 0, dims=2) .== 0) .* factorial.(R) .* prod(wbar .^ nubar, dims=2))
end

function probres(x, y, nu) # set input types
    N = 3
    out = Array{Float64, 1}(undef, length(x))
    w = [x y (1 .-x .-y)]
    flnw = floor.(N * w)
    nubar = nu .- flnw
    wbar = N * w .- flnw
    R = convert.(Int64, round.(sum(wbar,dims=2)))
    for i in 1:length(x)
        wbar[i,:] = (R[i]>0) ? (wbar[i,:] / R[i]) : (wbar[i,:] * 0)
        out[i] = (sum(nubar[i,:] .< 0) == 0) ? factorial(R[i]) * prod(wbar[i,:] .^ nubar[i,:]) : 0
    end
    return(out)
    #return((sum(nubar .< 0, dims=2) .== 0) .* factorial.(R) .* prod(wbar .^ nubar, dims=2))
end

using Plots

ternaryheatmap(prob111mn, 4)
ternaryheatmap(prob111res, 48)
