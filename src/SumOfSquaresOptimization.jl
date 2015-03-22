module SumOfSquaresOptimization

using SemidefiniteProgramming

# sos.jl
export SoS
export constraint, partconstraint, objective, partobjective
export @constraint, @partconstraint, @objective, @partobjective
export perturb_objective, perturb_objective_sparse

# solve.jl
export sossolve

# solution.jl
export SoSSolution, moment, @moment, dump, sosobj

# symmetrize.jl
export symmetrize!, symmetrize_cyclic!, symmetrize_dihedral!, symmetrize_full!

include("monoms.jl")
include("sos.jl")
include("solve.jl")
include("solution.jl")
include("symmetrize.jl")

end # module SumOfSquaresOptimization

