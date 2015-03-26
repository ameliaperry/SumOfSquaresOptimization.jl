module SumOfSquaresOptimization

importall SemidefiniteProgramming

# sos.jl
export Program
export constraint, partconstraint, objective, partobjective
export @constraint, @partconstraint, @objective, @partobjective
export perturb_objective, perturb_objective_sparse

# solve.jl
export sossolve

# solution.jl
export SoSSolution
export dump, obj, primalobj, dualobj
export moment, @moment

# symmetrize.jl
#export symmetrize!, symmetrize_cyclic!, symmetrize_dihedral!, symmetrize_full!

include("monoms.jl")
include("program.jl")
include("solve.jl")
include("solution.jl")
#include("symmetrize.jl")

end # module SumOfSquaresOptimization

