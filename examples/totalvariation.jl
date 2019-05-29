using SinkhornDistance, LinearAlgebra, Plots
using Distributions

dist1 = Normal(0.4,0.1)
dist2 = Normal(0.6, 0.1)

dist1_array = normalize(pdf.(dist1, 0:0.05:1), 1)
dist2_array = normalize(pdf.(dist2, 0:0.05:1), 1)

function tv_distmat(dim)
    result = ones(dim,dim)
    zero_el = zero(result[1,1])
    for i=1:dim
        result[i,i] = zero_el
    end
    result
end

true_tv(a, b) = 0.5 * norm(a .- b, 1)


distmat = tv_distmat(size(dist1_array,1))
plan = sinkhorn_plan(distmat, dist1_array, dist2_array; ϵ=1e-2, rounds=100)
cost = sum(plan .* distmat)

true_tv(dist1_array, dist2_array)
