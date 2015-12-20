using SumOfSquaresOptimization
using Base.Test

objtol = 1e-4

function equals_upto_perm(a,b)
    if length(a) != length(b)
        return false
    end
    mapreduce(p -> b[p] == a, Base.OrFun(), permutations(1:length(b)))
end
@test equals_upto_perm([1,2],[2,1])

# UNIT TESTS

# monomials & polynomials
a = Dict( :x1 => 3, :x2 => 2 ) 
b = Dict( :x1 => 1, :x3 => 5 )
c = Dict( :x1 => 4, :x2 => 2, :x3 => 5 )
@test a * b == c
@test SumOfSquaresOptimization.deg(a) == 5
p = Dict(a=>1.0, b=>2.0)
@test Dict(a=>1.0) + Dict(b=>2.0) == p
@test SumOfSquaresOptimization.deg(p) == 6
@test equals_upto_perm(SumOfSquaresOptimization.monoms([:x1,:x2],2), [ Dict(:x1=>2), Dict(:x1=>1, :x2=>1), Dict(:x2=>2) ])
@test SumOfSquaresOptimization.monoms(0) == [SumOfSquaresOptimization.SoSMonom()]
@test SumOfSquaresOptimization.monoms(3) == Array{SumOfSquaresOptimization.SoSMonom,1}()
@test equals_upto_perm(SumOfSquaresOptimization.decomp(Dict(:x1=>4),3), [ (Dict(:x1=>3),Dict(:x1=>1)), (Dict(:x1=>2),Dict(:x1=>2)), (Dict(:x1=>1),Dict(:x1=>3)) ])

# parsing
@test SumOfSquaresOptimization.parsepoly(Program(), parse("x1^3 * x2^2")) == Dict(a => 1.0)
@test SumOfSquaresOptimization.parsepoly(Program(), parse("x1^3 * x2^2 + 2.0 * x1 * x3^5")) == p



# END-TO-END TESTS

# unconstrained minimization: a quadratic bowl
prog = Program(minimize=true)
objective(prog, "x^2 + y^2 - 1")
sol = sossolve(prog,2)
@test abs(primalobj(sol) + 1) < objtol

# unconstrained, unbounded maximization
prog = Program(maximize=true)
objective(prog, "x^2 + y^2 - 1")
sol = sossolve(prog,2)
@test isinf(primalobj(sol))
@test !signbit(primalobj(sol))

# min bisection in the 6-cycle
prog = Program(minimize=true)
for i in 1:6
    for j in [-1,1]
        jp = (i + j + 2*6 - 1) % 6 + 1
        @partobjective(prog, "x%d * (1-x%d)", i, jp)
    end
    @constraint(prog, "x%d^2 == x%d", i, i)
    @partconstraint(prog, :bisect, "x%d", i)
end
@partconstraint(prog, :bisect, "-%d", 6/2)

sol = sossolve(prog,2)
@test abs(primalobj(sol) - 1.5) < objtol

# test dihedral symmetry on the same program
symmetrize_dihedral!(prog, [symbol("x$i") for i in 1:6])
sol = sossolve(prog,2)
@test abs(primalobj(sol) - 1.5) < objtol
sol = sossolve(prog,4)
@test abs(primalobj(sol) - 2.0) < objtol


