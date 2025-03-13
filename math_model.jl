using LinearAlgebra

r = 2.0
K = 120
P0 = [20 20 20 20 20]
N = sum(P0)
l = exp(r*(1-N/K))/3
m = [[0 0 l l l];
     [1 0 0 0 0];
     [0 1 0 0 0];
     [0 0 1 0 0];
     [0 0 0 1 0]]



P_size = []
for i=1:100
    N = sum(sum(P0))
    push!(P_size, N)
    l = .9exp(r*(1-N/K))/3
    if i == 100
        print(l)
    end
    m = [[0 0 l l l];
         [1 0 0 0 0];
         [0 1 0 0 0];
         [0 0 1 0 0];
         [0 0 0 1 0]]
    P0 = P0 * m
end
using Plots

P0/K

[1 1 1 2/3 1/3] * m
plot(P_size)