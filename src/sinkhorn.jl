function sinkhorn_plan(dist_mat, a, b; ϵ=1e-1, rounds=2)
    # figure out types later
    # do we want to allow for a version of this with batches over a and b?
    K = exp.(-dist_mat / ϵ)
    v = one.(b)
    u = a ./ (K * v)
    v = b ./ (K' * u)

    for iter=1:rounds
        u = a ./ (K * v)
        v = b ./ (K' * u)
    end
    u .* K .* v'
end

function softmin(mat, ϵ; dims=1)
    zbar = minimum(mat, dims=dims)
    (zbar .- sum(exp.(-(mat .- zbar)/ϵ), dims=dims))[:]
end

#function sinkhorn_plan_softmin(dist_mat, a, b; ϵ=1e-3, rounds=1)
    #hist_dim = size(dist_mat, 1)
    #g = ones(hist_dim)
    #f = ones(hist_dim)
    #S = dist_mat .- f .- g'
    #f .= softmin(S, ϵ, dims=2) .- f .+ ϵ*log.(a)
    #S .= dist_mat .- f .- g'
    #g .= softmin(S, ϵ, dims=1) .- g .+ ϵ*log.(b)
    #S .= dist_mat .- f .- g'
    #for iter=1:rounds
        #f .= softmin(S, ϵ, dims=2) .- f .+ ϵ*log.(a)
        #S .= dist_mat .- f .- g'
        #g .= softmin(S, ϵ, dims=1) .- g .+ ϵ*log.(b)
        #S .= dist_mat .- f .- g'
    #end
    #Diagonal(exp.(f ./ ϵ)) * exp.(-dist_mat / ϵ) * Diagonal(exp.(g ./ ϵ))
#end

function sinkhorn_plan_log(dist_mat, a, b; ϵ=1e-1, rounds=2)
    K = exp.(-dist_mat / ϵ)
    g = one.(a)
    f = ϵ*log.(a) - ϵ*log.(K * exp.(g / ϵ))
    g = ϵ*log.(b) - ϵ*log.(K' * exp.(f / ϵ))

    iters = 0

    for iter=1:rounds
        f = ϵ*( log.(a) - log.(K * exp.(g / ϵ)) )
        g = ϵ*( log.(b) - log.(K' * exp.(f / ϵ)) )
    end

    exp.(f ./ ϵ) .* K .* exp.(g ./ ϵ)'
end
