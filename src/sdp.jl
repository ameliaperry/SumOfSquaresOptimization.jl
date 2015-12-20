
# procedure for SDP, *in this order*:
# instantiate SDPSession(nconstraints, nmoments)
# call sdp_obj! a bunch of times in order to set objective
# call sdp_con! a bunch of times (in any order) to set constraints)
# call sdp_solve to return a solution


# solution structure & methods
type SDPSolution
    primalobj :: Float64
    dualobj :: Float64
    primalmatrix :: Array{Float64,2}
    dualmatrix :: Array{Float64,2}
    status :: Symbol # :Normal, :Infeasible, :Unbounded, :Warning, :Error
end

function sdp_sol_infeasible(maximize :: Bool)
    obj = maximize ? -Inf : Inf
    SDPSolution(obj, obj, zeros(0,0), zeros(0,0), :Infeasible)
end

function sdp_sol_unbounded(maximize :: Bool)
    obj = maximize ? Inf : -Inf
    SDPSolution(obj, obj, zeros(0,0), zeros(0,0), :Unbounded)
end




type SDPSession
    maximize :: Bool
    nmoments :: Int64
    initval :: SparseMatrixCSC{Float64,Int64}
    nconstraints :: Int64
    constraints_missing_objective :: Int64
    objective :: Array{Float64,1}
    offset :: Float64
    fname :: AbstractString
    io :: IO
end
SDPSession(maximize, nconstraints, nmoments, offset) = SDPSession(maximize, nmoments, spzeros(nmoments,nmoments), nconstraints, nconstraints, Float64[], offset, sdp_init(nconstraints,nmoments)...)

function sdp_init(nconstraints, nmoments)
    fname,io = mktemp()
    println(io, nconstraints)
    println(io, 1) # number of blocks
    println(io, nmoments)
    fname,io
end

function sdp_obj!(sess :: SDPSession, val :: Float64)
    push!(sess.objective, val)
    @printf(sess.io, "%.18e ", sess.maximize ? -val : val)
    if (sess.constraints_missing_objective -= 1) == 0
        println(sess.io)
    end
end

function sdp_con!(sess :: SDPSession, con :: Int64, i :: Int64, j :: Int64, val :: Float64)
    if con > sess.nconstraints || con < 0
        throw(ArgumentError("invalid constraint number $con ($(sess.nconstraints) expected"))
    end
    if sess.nconstraints == 0 
        sess.initval[i,j] = val
        sess.initval[j,i] = val
    end
    @printf(sess.io, "%d 1 %d %d %.18e\n", con, i, j, val)
end


retcode_to_status = [
    :Unbounded, #"Success, but the problem is primal infeasible.",
    :Infeasible, #"Success, but the problem is dual infeasible.",
    :Warning, #"Partial success: full accuracy was not achieved.",
    :Error, #"Failure: maximum iterations reached.",
    :Error, #"Failure: stuck at edge of primal infeasibility.",
    :Error, #"Failure: stuck at edge of dual infeasibility.",
    :Error, #"Failure: lack of progress.",
    :Error, #"Failure: X, Z, or O was singular.",
    :Error, #"Failure: detected NaN on Inf values."
]

function sdp_solve(sess :: SDPSession; call="csdp")
    close(sess.io)

    # First handle totally unconstrained SDPs -- these throw an error in CSDP.
    # The dual matrix is totally constrained to minus `initval` -- check if this is PSD.
    if sess.nconstraints == 0 
        eigmin = eigs(-sess.initval; nev=1, which=:SR)[1][1] # least eigenvalue
        if eigmin > -1e-9 # initial value is PSD, feasible
            return SDPSolution(sess.offset,sess.offset,full(-sess.initval),zeros(sess.nmoments,sess.nmoments),0)
        else # infeasible
            println("The SDP is infeasible.")
            return sdp_unbounded(sess.maximize)
        end
    end

    outfname, outio = mktemp()

    primalobj = NaN
    dualobj = NaN
    retcode = 0

    # note: our primal is CSDP's dual, so everything here is flipped that way
    for l in eachline(`$call $(sess.fname) $outfname`)
        print(l)

        if startswith(l, "Failure: return code is ")
            retcode = int(strip(split(l, "is ")[2]))
            if retcode > 0 && retcode <= length(csdp_messages)
                println(csdp_messages[retcode])
            else
                println("Unrecognized return value.")
            end

        elseif startswith(l, "Success: SDP is primal infeasible")
            println("The SDP is infeasible.")
            return sdp_sol_infeasible(sess.maximize)
        elseif startswith(l, "Success: SDP is dual infeasible")
            println("The SDP is unbounded.")
            return sdp_sol_unbounded(sess.maximize)
        end
    end

    dualvector = map(float, split(chomp(readline(outio))))

    primalmatrix = zeros(sess.nmoments, sess.nmoments)
    dualmatrix = zeros(sess.nmoments, sess.nmoments)
    for l in eachline(outio)
        toks = split(chomp(l), " ")
        i = parse(Int,toks[3])
        j = parse(Int,toks[4])
        v = parse(Float64,toks[5])
        if parse(Int,toks[1]) == 1
            primalmatrix[i,j] = v
            primalmatrix[j,i] = v
        else
            dualmatrix[i,j] = v
            dualmatrix[j,i] = v
        end
    end
    
    # compute objective
    if retcode != 1 && retcode != 2
        primalobj = dot(dualvector, sess.objective) + sess.offset
        dualobj = primalobj # XXX this is no good -- want to know the duality gap.
    end

    SDPSolution(primalobj, dualobj, primalmatrix, dualmatrix, retcode==0 ? :Normal : retcode_to_status[retcode])
end

