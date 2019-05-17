using SinkhornDistance, LinearAlgebra, Plots
using Distributions

function example_data(n)
    distmat = zeros(n,n)
    marginal_a = [exp(-(x - n/2)^2/8) for x=1:n]
    marginal_b = [exp(-(x - 3n/4)^2/8) for x=1:n]
    marginal_b .+= [exp(-(x - n/4)^2/8) for x=1:n]
    marginal_a ./= sum(marginal_a)
    marginal_b ./= sum(marginal_b)
    for i=1:n
        for j=1:n
            distmat[i, j] = abs((i/n)-(j/n))
        end
    end
    (marginal_a, marginal_b, distmat)
end

a, b, distmat = example_data(30)
standard_plan = sinkhorn_plan(distmat, a, b; ϵ=1e-2, rounds=4)

log_plan = sinkhorn_plan_log(distmat, a, b; ϵ=1e-2)

s, _ = qr(rand(10,10))
