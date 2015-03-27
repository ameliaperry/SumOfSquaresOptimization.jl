# SumOfSquaresOptimization

A Julia library to solve sum-of-squares relaxations of polynomial systems. This is not yet in the Julia package repository.

## example

The following solves the degree 4 SoS relaxation for vertex cover on the complete graph on five vertices:

``` julia
using SumOfSquaresOptimization

prog = Program(minimize=true)

for i in 1:5
    for j in 1:5
        if j > i
            @constraint(prog, "x%d + x%d >= 1", i,j) # coverage constraint
        end
    end
    @constraint(prog, "x%d^2 - x%d = 1", i,i) # hypercube constraint
    @partobjective(prog, "x%d", i) # objective: min sum x_i
end

sol = sossolve(prog,4) # 4 indicates degree
dump(sol)
```

## functions

The following four macros specify constraints and objectives in `printf` style:
* `@constraint(prog, fmt, ...)` adds a polynomial constraint
* `@objective(prog, fmt, ...)` sets a polynomial objective
* `@partconstraint(sos, key, fmt, ...)` adds a polynomial to the constraint labeled by `key`
* `@partobjective(sos, fmt, ...)` adds a polynomial to the current objective
There are also function analogues for plain strings: `constraint`, `objective`, `partconstraint`, and `partobjective`.

The following will solve a polynomial system and access a solution:
* `sossolve(sos, deg)` solves the program, returning a solution object. You can optionally specify `solver="csdp"` or `solver="sdpa"`, with `csdp` being the default. One of these two must be installed.
* `dump(sol)` outputs the objective, moments, and dual matrix to `stdout`. These will be more readable in the future.
* `@moment(sol, fmt, ...)` outputs the pseudo-expectation of the given polynomial, specified in `printf` style.
* `primalobj(sol)` and `dualobj(sol)` return the objective values (which are hopefully equal).

The following functions provide symmetry hints, enabling faster solution and ensuring symmetric moments:
* `symmetrize!(prog, perms)`, where `perms` is an array of `Dict{Symbol,Symbol}` objects encoding *generating* permutations.
* `symmetrize_dihedral!(prog, cycle)`, where `cycle` is an array of `Symbol`s. This imposes symmetry under the action of the dihedral group.
* `symmetrize_cyclic!(prog, cycle)`, much as the previous but for the cyclic group -- so no reflection symmetry of the given cycle is imposed.
* `symmetrize_full!(prog, symbols)`, for full permutation symmetry of the given `Array{Symbol}`.


## troubleshooting
There are probably still a lot of issues with this code, and feel free to let me know about them. This package is in very early stages still.

One issue I've encountered is that `csdp` doesn't seem to accept more than 23169 constraints in some compiled 32-bit versions, and these SoS programs quickly excced that (e.g. by stepping up the degree to 6 in the example above).

## todo
* actually test `sdpa` support
* pretty-print the moments
* pretty-print Psatz refutation
* write tests
* investigate the JuliaOpt SDP packages and compare
* any way to do high-precision solution, leveraging Julia's `BigFloat` and multi-precision SDPA?

