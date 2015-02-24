# SumOfSquaresOptimization

[![Build Status](https://travis-ci.org/willperry/SumOfSquaresOptimization.jl.svg?branch=master)](https://travis-ci.org/willperry/SumOfSquaresOptimization.jl)

A Julia library to solve sums-of-squares relaxations of polynomial systems.

## example

The following solves the degree 2 SoS relaxation for vertex cover on the complete graph on five vertices:

``` julia
using SumOfSquaresOptimization

sos = SoS(minimize=true)

for i in 1:5
    for j in 1:5
        if j > i
            @constraint(sos, "x%d + x%d >= 1", i,j) # coverage constraint
        end
    end
    @constraint(sos, "x%d^2 - x%d == 1", i,i) # hypercube constraint
    @partobjective(sos, "x%d", i) # objective: min sum x_i
end

sol = sossolve(sos,2) # 2 indicates degree
```



