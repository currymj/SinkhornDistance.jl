using Flux
using SinkhornDistance
using Plots

fmodel = Chain(Dense(10,128),
Dense(128,128),
Dense(128,100),
x->reshape(x, 10,10).^2,
x->x'*x)

a = rand(10)
a ./= sum(a)
b = rand(10)
b ./= sum(b)
distmat = fmodel(rand(10))

diffplan = sinkhorn_plan(distmat, a, b)

resultloss = sum(distmat .* diffplan)
