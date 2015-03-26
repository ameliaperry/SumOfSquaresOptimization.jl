
# The output of the SoS method on a polynomial program

type SoSSolution
    prog :: Program
    degree :: Int64
    primalobj :: Float64
    dualobj :: Float64
    moments :: Dict{SoSMonom,Float64}
    dualmatrix :: SemidefiniteProgramming.SparseSymmetricBlockMatrix{Float64}
end

function obj(sol :: SoSSolution)
    sol.primalobj
end
function primalobj(sol :: SoSSolution)
    sol.primalobj
end
function dualobj(sol :: SoSSolution)
    sol.dualobj
end


# Extracts moments from a solution.

function moment(sol :: SoSSolution, monom :: SoSMonom)

    if deg(monom) > sol.degree
       throw(ArgumentError("Requested moment degree is higher than solution degree"))
    elseif !haskey(sol.moments, monom)
       throw(ArgumentError("Requested moment is not described by this SoS solution"))
    end

    sol.moments[monom]
end

function moment(sol :: SoSSolution, poly :: SoSPoly)
    tot = 0.0
    for (monom,coeff) in poly
        tot += moment(sol, monom) * coeff
    end
    tot
end

function moment(sol :: SoSSolution, poly :: Expr)
    moment(sol, parsepoly(nothing, poly))
end

function moment(sol :: SoSSolution, symb :: Symbol)
    moment(sol, parsepoly(nothing, symb))
end

macro moment(sol, args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        moment( $(esc(sol)), parse(str) )
    end
end


# TODO: some sort of handle on the dual solution.


# Dumps a solution to stdout.

function dump(sol :: SoSSolution)
    println("Objective:")
    println(sol.primalobj)
    println("Moments:")
    println(sol.moments)
    println("Dual matrix:")
    println(sol.dualmatrix)
end

