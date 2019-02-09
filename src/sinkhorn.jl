function sinkhorn_plan(dist_mat, a, b; ϵ=1e-3, tol=1e-3)
    # figure out types later
    # do we want to allow for a version of this with batches over a and b?
    K = exp.(-dist_mat / ϵ)
    K ./= sum(K)
    hist_dim = size(K, 1)
    v = ones(hist_dim)
    u = a ./ (K * v)
    v .= b ./ (K * u)
    prev_u = ones(hist_dim)

    iters = 0
    # michiel stock's stopping criterion
    while maximum(abs.(prev_u - u)) >= tol# switch this to test convergence instead
        prev_u .= u
        u .= a ./ (K * v)
        v .= b ./ (K * u)
        iters = iters + 1
    end
    println(iters)
    Diagonal(u) * K * Diagonal(v)
end
