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
dump(sol)
```

## functions

* `@constraint(sos, fmt, ...)` adds a polynomial constraint, specified in `printf` style.
* `@objective(sos, fmt, ...)` sets a polynomial objective, specified in `printf` style.
* `@partobjective(sos, fmt, ...)` adds a polynomial to the current objective, specified in `printf` style.
* `sossolve(sos, deg)` solves the program, returning a solution object. You can optionally specify `solver="csdp"` or `solver="sdpa"`, with `csdp` being the default. One of these two must be installed.
* `dump(sol)` outputs the objective, primal matrix, and dual matrix. These will be more readable in future releases.
* `@moment(sol, fmt, ...)` outputs the pseudo-expectation of the given polynomial, specified in `printf` style.

## troubleshooting
There are probably still a lot of issues with this code, and feel free to let me know about them. This package is in very early stages still.

One issue I've encountered is that `csdp` doesn't seem to accept more than 23169 constraints in some compiled 32-bit versions, and these SoS programs quickly excced that (e.g. by stepping up the degree to 6 in the example above).

## todo
* actually test `sdpa` support
* pretty-print the moments
* pretty-print Psatz refutation
* write tests
* investigate the JuliaOpt SDP packages and compare
* investigate ways to produce a less huge SDP

