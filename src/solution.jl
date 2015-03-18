

type SoSSolution
    prog :: SoS
    degree :: Int64
    alphaidx :: Dict{(SoSMonom,SoSMonom,SoSMonom,SoSMonom),Int64}
    betaidx :: Dict{(SoSPoly,SoSMonom),Int64}
    sdp :: SparseSDPSolution
end

function sosobj(sol :: SoSSolution)
    obj(sol.sdp)
end

function moment(sol :: SoSSolution, monom :: SoSMonom)

    if(deg(monom) > sol.degree)
        throw(ArgumentError("Requested moment degree is higher than solution degree"))
    end

    (m1,m2) = decomp1(monom, div(sol.degree,2))
    
    mat = primalmatrix(sol.sdp)[1]
    # TODO we should complain if the variables in m1 or m2 aren't variables of our SoS program.
    # Currently unrecognized moments return 0 (this is the implementation in SemidefiniteProgramming.jl)
#    if haskey(entries(mat), (m1,m2))
#        throw(ArgumentError("Solution does not have moments for these variables"))
#    end
    mat[m1,m2]
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

function dump(sol :: SoSSolution)
    println("Objective:")
    println(obj(sol.sdp))
    println("Primal:")
    println(primalmatrix(sol.sdp))
    println("Dual:")
    println(dualmatrix(sol.sdp))
end

