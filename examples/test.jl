using SinkhornDistance, LinearAlgebra, Plots

a = ones(10)
a ./= norm(a)
b = ones(10)
b ./= norm(b)
distmat = rand(10,10)
distmat = distmat' * distmat
distmat ./= norm(distmat)
standard_plan = sinkhorn_plan(distmat, a, b; rounds=3)

log_plan = sinkhorn_plan_log(distmat, a, b; ϵ=1e-2)
