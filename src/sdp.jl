

# procedure for SDP, *in this order*:
# instantiate SDPSession(nconstraints, nmoments)
# call sdp_obj! a bunch of times in order to set objective
# call sdp_con! a bunch of times (in any order) to set constraints)
# call sdp_solve to return a solution

type SDPSession
    maximize :: Bool
    nmoments :: Int64
    nconstraints :: Int64
    constraints_missing_objective :: Int64
    objective :: Array{Float64,1}
    offset :: Float64
    fname :: String
    io :: IO
end
SDPSession(maximize, nconstraints, nmoments, offset) = SDPSession(maximize, nmoments, nconstraints, nconstraints, Float64[], offset, sdp_init(nconstraints,nmoments)...)

type SDPSolution
    primalobj :: Float64
    dualobj :: Float64
    primalmatrix :: Array{Float64,2}
    dualmatrix :: Array{Float64,2}
    retcode :: Int
end

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
    @printf(sess.io, "%d 1 %d %d %.18e\n", con, i, j, val)
end


csdp_messages = [
    "Success, but the problem is primal infeasible.",
    "Success, but the problem is dual infeasible.",
    "Partial success: full accuracy was not achieved.",
    "Failure: maximum iterations reached.",
    "Failure: stuck at edge of primal infeasibility.",
    "Failure: stuck at edge of dual infeasibility.",
    "Failure: lack of progress.",
    "Failure: X, Z, or O was singular.",
    "Failure: detected NaN on Inf values."
]

function sdp_solve(sess :: SDPSession)
    close(sess.io)
    outfname, outio = mktemp()

    primalobj = NaN
    dualobj = NaN
    retcode = 0

    # note: our primal is CSDP's dual, so everything here is flipped that way
    for l in eachline(`csdp $(sess.fname) $outfname`)

        if beginswith(l, "Failure: return code is ")
            retcode = int(strip(split(l, "is ")[2]))
            if retcode > 0 && retcode <= length(csdp_messages)
                println(csdp_messages[retcode])
            else
                println("Unrecognized return value.")
            end

        elseif beginswith(l, "Success: SDP is primal infeasible")
            primalobj = sess.maximize ? Inf : -Inf
            dualobj = sess.maximize ? Inf : -Inf
            retcode = 1
            println("The program is dual infeasible.")
        
        elseif beginswith(l, "Success: SDP is dual infeasible")
            primalobj = sess.maximize ? -Inf : Inf
            dualobj = sess.maximize ? -Inf : Inf
            retcode = 2
            println("The program is primal infeasible.")

        end
    end

    # compute objective
    if retcode != 1 && retcode != 2

    end
    
    dualvector = map(float, split(chomp(readline(outio))))
    primalobj = dot(dualvector, sess.objective) + sess.offset

    primalmatrix = zeros(sess.nmoments, sess.nmoments)
    dualmatrix = zeros(sess.nmoments, sess.nmoments)
    for l in eachline(outio)
        toks = split(chomp(l), " ")
        i = int(toks[3])
        j = int(toks[4])
        v = float(toks[5])
        if int(toks[1]) == 1
            primalmatrix[i,j] = v
        else
            dualmatrix[i,j] = v
            dualmatrix[j,i] = v
        end
    end
    dualobj = primalobj # XXX this is no good. should we take the truncated CSDP output?

    SDPSolution(primalobj, dualobj, primalmatrix, dualmatrix, retcode)
end

