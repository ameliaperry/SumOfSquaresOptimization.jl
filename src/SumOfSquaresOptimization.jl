module SumOfSquaresOptimization
import Base.Operators: +, *
using DataStructures # for SortedDict, for row-redundancy checks
import DataStructures: status # annoying that we have to overload this with a different meaning

# sos.jl
export Program
export constraint, partconstraint, objective, partobjective
export @constraint, @partconstraint, @objective, @partobjective
export perturb_objective, perturb_objective_sparse

# solve.jl
export sossolve

# solution.jl
export SoSSolution
export dumpsol, obj, primalobj, dualobj, status
export moment, @moment

# symmetrize.jl
export symmetrize!
export symmetrize_cyclic!, symmetrize_dihedral!, symmetrize_hypercube!, symmetrize_full!

include("monoms.jl")
include("program.jl")
include("sdp.jl")
include("solve.jl")
include("solution.jl")
include("symmetrize.jl")

end # module SumOfSquaresOptimization

