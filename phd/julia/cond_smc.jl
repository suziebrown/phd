function csmc(N::Int64, T::Int64, observations::Array{Float64,1} initialsam::Function, transition::Function, potential::Function, immortal_positions::Array{Float64,1})
    # pre-allocate memory
    immortal_indices = Array{Int64, 1}(undef, T+1)
    positions_old = Array{Float64, 1}(undef, N) # may not be necessary to store these
    positions_new = Array{Float64, 1}(undef, N) # probably doesn't need initialising
    # (also need to store ancestry tree but haven't worked out an efficient way to update it yet.)

    # set OU process parameters
    delta = 1.0
    sigma = 0.1

    # initialise
    positions_old = initialsam(N, delta, sigma)
    immortal_indices = rand(1:N, T+1) # pre-sample all of them to speed up computation
end
