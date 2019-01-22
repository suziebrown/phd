using Random, Distributions

function ousim(T::Int64, delta::Float64, sigma::Float64, returnstates::Bool)
    # pre-allocate arrays
    states = Array{Float64}(undef, T+1)
    observations = Array{Float64}(undef, T+1)

    # generate Normal noise:
    states_noise = rand(Normal(0, delta^(0.5)), T+1)
    obs_noise = rand(Normal(0, sigma), T+1)

    # initial states
    states[1] = states_noise[1]
    observations[1] = states[1] + obs_noise[1]

    # simulate the process
    for t in 2:T+1
        states[t] = (1-delta) * states[t-1] + states_noise[t]
        observations[t] = states[t] + obs_noise[t]
    end

    if returnstates
        return (states = states, observations = observations)
    else
        return observations
    end

end


# TEST
using Plots
gr()

x = ousim(10, 0.1, 1.0, true)
plot(0:10, x.states, label="states")
plot!(0:10, x.observations, label="observations")

x = ousim(10, 0.1, 1.0, false)
plot(0:10, x, label="observations")
