module SumOfSquaresOptimization

using SemidefiniteProgramming

export SoS, constraint, partconstraint, objective, partobjective, @constraint, @partconstraint, @objective, @partobjective
export sossolve
export SoSSolution, moment, @moment, dump, sosobj
export symmetrize!, symmetrize_cyclic!, symmetrize_dihedral!, symmetrize_full!

include("monoms.jl")
include("sos.jl")
include("solve.jl")
include("solution.jl")
include("symmetrize.jl")

end # module SumOfSquaresOptimization

