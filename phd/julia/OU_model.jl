using Random, Distributions

# generate a sequence of observations from an Orstein-Uhlenbeck process
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

# sample initial particle positions
function ouinit(N::UInt16, delta::Float64, sigma::Float64)
    rand(Normal(0, delta^(0.5)), N)
end

# propagate states: sample new positions from current ones
function outransition(oldpos::Array{Float64,1}, delta::Float64, sigma::Float64)
    N = length(oldpos)
    rand(Normal(), N) .* delta^(0.5) .+ (1-delta) .* oldpos
end

# calculate potentials between particle positions and observations (to compute weights)
function oupotential(pos::Array{Float64,1}, obs::Float64, delta::Float64, sigma::Float64)
    (2*pi)^(-0.5) * sigma^(-1) * exp.(-(obs .- pos).^2 / (2*sigma^2))
end
