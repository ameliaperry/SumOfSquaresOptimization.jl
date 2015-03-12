

type SoSSolution
    degree :: Int64
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
    
    # XXX what if the input variables are bad, and we don't have moments for them?
    # should test that primalmatrix actually has these entries, and throw a 
    # meaningful exception if not
    primalmatrix(sol.sdp)[1][m1,m2]
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

