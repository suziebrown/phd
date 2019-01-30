# kalman filter & rts smoother

function oukalman(delta::Float64, sigma::Float64, observations::Array{Float64,1})
    T = length(observations) - 1

    # assign memory
    xhat = Array{Float64, 1}(undef, T+1) # means \hat{x}_{k|k}
    sigma = Array{Float64,1}(undef, T+1) # variances \Sigma_{k|k}

    # initialise
    tempx = 0
    tempsigma = 1

    # recursion
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

function ourts(delta::Float64, sigma::Float64, observations::Array{Float64,1})
    # call to oukalman in here
    #NOT FINISHED
end
