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
        # NOT FINISHED
    end
end

function ourts(delta::Float64, sigma::Float64, observations::Array{Float64,1})
    # call to oukalman in here
    #NOT FINISHED
end
