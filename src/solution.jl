
# The output of the SoS method on a polynomial program

type SoSSolution
    prog :: Program
    degree :: Int64
    primalobj :: Float64
    dualobj :: Float64
    moments :: Dict{SoSMonom,Float64}
    primalmatrix :: Array{Float64,2}
    dualmatrix :: Array{Float64,2}
    status :: Symbol # :Normal, :Infeasible, :Unbounded, :Warning, :Error
end

obj(sol :: SoSSolution) = sol.primalobj
primalobj(sol :: SoSSolution) = sol.primalobj
dualobj(sol :: SoSSolution) = sol.dualobj
status(sol :: SoSSolution) = sol.status

function sos_sol_infeasible(prog :: Program, deg :: Int)
    obj = prog.maximize ? -Inf : Inf
    SoSSolution(prog,deg,obj,obj,Dict{SoSMonom,Float64}(),zeros(0,0),zeros(0,0),:Infeasible)
end
function sos_sol_unbounded(prog :: Program, deg :: Int)
    obj = prog.maximize ? Inf : -Inf
    SoSSolution(prog,deg,obj,obj,Dict{SoSMonom,Float64}(),zeros(0,0),zeros(0,0),:Unbounded)
end



function badsol_check(sol :: SoSSolution)
    if sol.status == :Unbounded
        throw(ArgumentError("The program was unbounded."))
    elseif sol.status == :Infeasible
        throw(ArgumentError("The program was infeasible."))
    elseif sol.status == :Error
        throw(ArgumentError("There was an error in solving the program."))
    elseif sol.status == :Warning
        println("Warning: the SDP solution was abnormal in some way.")
    end
end


# Extracts moments from a solution.

function moment(sol :: SoSSolution, monom :: SoSMonom)

    if deg(monom) > sol.degree
       throw(ArgumentError("Requested moment degree is higher than solution degree"))
    elseif !haskey(sol.moments, monom)
       throw(ArgumentError("Requested moment is not described by this SoS solution"))
    end
    badsol_check(sol)

    sol.moments[monom]
end

moment(sol :: SoSSolution, poly :: SoSPoly) = sum([ moment(sol,monom) * coeff for (monom,coeff) in poly ])
moment(sol :: SoSSolution, poly :: Expr) = moment(sol, parsepoly(nothing, poly))
moment(sol :: SoSSolution, symb :: Symbol) = moment(sol, parsepoly(nothing, symb))
moment(sol :: SoSSolution, str :: AbstractString) = moment(sol, parse(str))
macro moment(sol, args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        moment( $(esc(sol)), parse(str) )
    end
end




# Polynomial represented by the dual matrix. This is congruent,
# modulo the constraints, to the objective function minus the dual
# objective value, and is a sum of squares.
# XXX currently unused and private
function dualpoly(sol :: SoSSolution)
    mon = vcat( [ monoms(sol.prog,i) for i in 0:div(sol.degree,2) ] ... )

    po = SoSPoly()
    for i in 1:length(mon)
        for j in 1:length(mon)
            prod = mon[i] * mon[j]
            po[prod] = get(po, prod, 0.0) + sol.dualmatrix[i,j]
        end
    end
    po
end


# This code processes the dual solution / SoS proof, and prints polynomials, the sum of
# whose squares equals the dual polynomial (cf dualpoly()).
# XXX currently unused, private, and broken
function printsquares(sol :: SoSSolution; adjust=1.0, method=:svd, thresh=0.00001)

    # decode the dual matrix
    dm = sol.dualmatrix
    #mon = [ SumOfSquaresOptimization.monoms(prog,i) for i in 0:div(deg,2) ]
    mon = [ SumOfSquaresOptimization.monoms(prog,i) for i in div(deg,2):-1:0 ] # reverse order for good cholesky results
    mon = vcat(mon...)
    mtx = zeros(length(mon),length(mon))
    for i in 1:length(mon)
        for j in 1:length(mon)
            mtx[i,j] = dm[1,mon[i],mon[j]]
        end
    end

    if method == :cholesky
        ch=chol(mtx)

        nonzero = 0
        for i in 1:size(mtx,1)
            vec = Base.vec(ch[i,:])

            # threshold
            if dot(vec,vec) < thresh
                continue
            end

            # make the poly
            s = SoSPoly()
            for j in 1:size(mtx,1)
                if abs(vec[j]) < thresh
                    continue
                end
                s[mon[j]] = vec[j] / sqrt(adjust)
            end

            # display
            @printf("square of: %s\n\n", s)
            nonzero += 1

        end
        @printf("threshold after %d squares; the remaining %d are negligible\n", nonzero, size(mtx,1)-nonzero)

    else
        # svd to get the squares involved
        sv=svdfact(mtx)

        for i in 1:size(mtx,1)
            # singular value
            sc = sqrt(sv[:S][i] / adjust)

            # threshold
            if sc^2 < thresh
                @printf("threshold after %d squares; the remaining %d are negligible\n", i-1, size(mtx,1)-i+1)
                break
            end

            # make the poly
            s = SoSPoly()
            for j in 1:size(mtx,1)
                coeff = sv[:U][j,i] * sc
                if(abs(coeff) >= thresh)
                    s[mon[j]] = coeff
                end
            end

            # display
            @printf("square of: %s\n\n", s)
        end
    end
end


# Dumps a solution to stdout.
function dumpsol(sol :: SoSSolution)
    badsol_check(sol)
    println("Objective:")
    println(sol.primalobj)
    println("Moments:")
    println(sol.moments)
    println("Dual matrix:")
    println(sol.dualmatrix)
end

