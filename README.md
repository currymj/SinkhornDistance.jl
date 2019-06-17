# SinkhornDistance.jl

a work in progress of the Sinkhorn matrix scaling algorithm for regularized
optimal transport. The plan is to have a small, simple implementation that
works on both CPU and GPU, and with autodiff systems. currently, it seems to work for this purpose (compatible with CuArrays and Flux). There is likely room for optimization, however.

The functions `sinkhorn_plan` and `sinkhorn_plan_log` take in a distance matrix and two marginals and spit out the optimal regularized plan after some number of iterations.

Possible future improvements include:

- convenience methods to directly compute the Sinkhorn distance as a single number
- convenience methods using [Distances.jl](https://github.com/JuliaStats/Distances.jl) to create distance matrices for common ground costs (Euclidean/Wasserstein, 0-1/Total Variation, etc.)

# relevant papers

The implementations are taken from Ch. 4 of the excellent:

Gabriel Peyré and Marco Cuturi. 2018. Computational Optimal Transport. arXiv:1803.00567 \[stat\] (March 2018). Retrieved October 31, 2018 from http://arxiv.org/abs/1803.00567

The paper that introduced the use of the Sinkhorn algorithm for optimal transport is:

Marco Cuturi. 2013. Sinkhorn distances: Lightspeed computation of optimal transport. In Advances in neural information processing systems, 2292–2300.


