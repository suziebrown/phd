using Random, Distributions

function ousim(T, delta, sigma)
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

    return (states = states, observations = observations)
end


# TEST
using Plots
gr()

x = ousim(10, 0.1, 1)
plot(0:10, x.states, label="states")
plot!(0:10, x.observations, label="observations")
