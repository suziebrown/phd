# kalman filter & rts smoother

function oukalman(delta::Float64, sigma::Float64, observations::Array{Float64,1})
    T = length(observations) - 1

    # assign memory PREDICTS CAN BE OVERWRITTEN EACH ITERATION
    predict = Array{Float64, 1}(undef, T+1) # x_{k|k-1}
    filter = Array{Float64, 1}(undef, T+1) # x_{k|k}
    predictvar = Array{Float64, 1}(undef, T+1) # sigma_{k|k-1}
    filtervar = Array{Float64,1}(undef, T+1) # sigma_{k|k}

    # initialise
    predict[1] = 0
    predictvar[1] = 1

    # recursion
    for t in 1:T
        a = predictvar[t] / (predictvar[t] + sigma^2)
        filter[t] = predict[t] + a * (observations[t] - predict[t])
        filtervar[t] = predictvar[t] * (1 - a)
        predict[t+1] = (1-delta) * filter[t]
        predictvar[t+1] = (1-delta)^2 * filtervar[t] + delta
    end

    # final filter state
    a = predictvar[T+1] / (predictvar[T+1] + sigma^2)
    filter[T+1] = predict[T+1] + a * (observations[T+1] - predict[T+1])
    filtervar[T+1] = predictvar[T+1] * (1 - a)
end

function ourts(delta::Float64, sigma::Float64, observations::Array{Float64,1})
    # call to oukalman in here
    #NOT FINISHED
end
