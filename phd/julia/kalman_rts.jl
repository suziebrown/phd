# Kalman filter & RTS smoother for OU model
#

function oukalman(delta::Float64, sigmain::Float64, observations::Array{Float64,1})
    T = length(observations) - 1

    # assign memory
    xhat = Array{Float64, 1}(undef, T+1) # means \hat{x}_{k|k}
    sigma = Array{Float64,1}(undef, T+1) # variances \Sigma_{k|k}

    # initialise
    tempx = 0 # \hat{x}_{k+1|k}
    tempsigma = 1 # \Sigma_{k+1|k}

    # forward recursion
    for t in 1:T
        a = tempsigma / (tempsigma + sigmain^2)
        xhat[t] = tempx + a * (observations[t] - tempx)
        sigma[t] = tempsigma * (1 - a)
        tempx = (1-delta) * xhat[t]
        tempsigma = (1-delta)^2 * sigma[t] + delta
    end

    # final filter state
    a = tempsigma / (tempsigma + sigmain^2)
    xhat[T+1] = tempx + a * (observations[T+1] - tempx)
    sigma[T+1] = tempsigma * (1 - a)

    return (mean = xhat, variance = sigma)
end

# could be useful to have a different version that takes kalman output as input...
function ourts(delta::Float64, sigmain::Float64, observations::Array{Float64,1})
    # forward Kalman pass
    kalman = oukalman(delta, sigmain, observations)
    xhat = kalman.mean
    sigma = kalman.variance
    # Warning: this function will incrementally OVERWRITE this Kalman output
    # (xhat, sigma) with the RTS output

    # backward recursion
    for i in 1:T
        t = T+1 - i # reverse time
        a = (1-delta) * sigma[t] / ((1-delta)^2 * sigma[t] + delta)
        xhat[t] = a * (xhat[t+1] - (1-delta) * xhat[t])
        sigma[t] = a * (a * sigma[t+1] - (1-delta) * sigma[t])
    end

    return (mean = xhat, variance = sigma)
end
