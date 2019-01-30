# Kalman filter & RTS smoother for OU model
#

function oukalman(delta::Float64, sigma::Float64, observations::Array{Float64,1})
    T = length(observations) - 1

    # assign memory
    xhat = Array{Float64, 1}(undef, T+1) # means \hat{x}_{k|k}
    sigma = Array{Float64,1}(undef, T+1) # variances \Sigma_{k|k}

    # initialise
    tempx = 0
    tempsigma = 1

    # forward recursion
    for t in 1:T
        a = tempsigma / (tempsigma + sigma^2)
        xhat[t] = tempx + a * (observations[t] - tempx)
        sigma[t] = tempsigma * (1 - a)
        tempx = (1-delta) * xhat[t]
        tempsigma = (1-delta)^2 * sigma[t] + delta
    end

    # final filter state
    a = tempsigma / (tempsigma + sigma^2)
    xhat[T+1] = tempx + a * (observations[T+1] - tempx)
    sigma[T+1] = tempsigma * (1 - a)

    return (mean = xhat, variance = sigma)
end

# could be useful to have a different version that takes kalman output as input...
function ourts(delta::Float64, sigma::Float64, observations::Array{Float64,1})
    # forward Kalman pass
    kalman = oukalman(delta, sigma, observations)
    kalmean = kalman.mean
    kalvar = kalman.variance

    # assign memory - ACTUALLY, we could overwrite kalmean/kalvar instead of saving both...
    xhat = Array{Float64, 1}(undef, T+1) # means \hat{x}_k
    sigma = Array{Float64,1}(undef, T+1) # variances \Sigma_k

    # initialise
    xhat[T+1] = kalmean[T+1]
    sigma[T+1] = kalvar[T+1]

    # backward recursion
    for i in 1:T
        t = T+1 - i # go backwards in time
        a = (1-delta) * kalvar[t] / ((1-delta)^2 * kalvar[t] + delta)
        xhat[t] = a * (xhat[t+1] - (1-delta) * kalmean[t])
        sigma[t] = a * (a * sigma[t+1] - (1-delta) * kalvar[t])
    end

    return (mean = xhat, variance = sigma)
end
