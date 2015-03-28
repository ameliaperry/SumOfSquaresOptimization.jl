using SumOfSquaresOptimization
using Base.Test

objtol = 1e-4

# UNIT TESTS

# TODO

# END-TO-END TESTS

# unconstrained minimization
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
symmetrize_dihedral!(prog, [symbol("x$i") for i in 1:6])
sol = sossolve(prog,2)
@test abs(primalobj(sol) - 1.5) < objtol
sol = sossolve(prog,4)
@test abs(primalobj(sol) - 2.0) < objtol


