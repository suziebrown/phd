N = 100
nrep = 100
k = Array{Int64,1}(undef, nrep)

for i in 1:nrep
    w = rand(N)
    w = w/sum(w)
    k[i] = N - sum(floor.(N*w))
end

histogram(k, title="N=100")

meank = Array{Float64,1}(undef, 1000)
vark = Array{Float64,1}(undef, 1000)
for N in 1:1000
    for i in 1:nrep
        w = rand(Normal(),N)
        w = w/sum(w)
        k[i] = N - sum(floor.(N*w))
    end
    meank[N] = mean(k)
    vark[N] = var(k)
end

plot(1:1000, meank, ribbon=(vark.^0.5, vark.^0.5), fillalpha=0.3, leg=false, title="Normal weights", xlabel="N", ylabel="k")

plot(1:1000, vark)
