using SinkhornDistance, Test, LinearAlgebra

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
function standard_domain_test()
    a, b, distmat = example_data(20)
    sinkhorn_plan(distmat, a, b; rounds=3)
end

function log_domain_test()
    a, b, distmat = example_data(20)
    sinkhorn_plan_log(distmat, a, b; rounds=1)
end

function tests()
    @testset "numerical stability on reasonable inputs" begin
        @test !any(isnan.(standard_domain_test()))
        @test !any(isnan.(log_domain_test()))
    end
end

tests()
