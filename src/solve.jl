#####################
###### Solving ######
#####################

#  Make the magic happen.
function sossolve(sos :: SoS, d :: Int64; solver="csdp")

    # sanity check
    if(d < 0)
        throw(ArgumentError("Degree must be non-negative"))
    elseif(d < deg(sos.objective))
        throw(ArgumentError("Degree must be at least that of the objective"))
    elseif(d % 2 == 1)
        throw(ArgumentError("Degree must be even"))
    end

    
    # set up problem & solver instance
    sdp = SparseSDP(maximize=sos.maximize)
    if (solver == "sdpa")
        solverinst = SDPA()
    elseif (solver == "csdp")
        solverinst = CSDP()
    else
        throw(ArgumentError("Unsupported solver."))
    end


    
    ## objective ##
    for (k,v) in sos.objective
        (k1,k2) = decomp1(k,div(d,2)) # just one representative for each monomial
        #           blk row col val
        setobj!(sdp, 1, k1, k2, v)
    end


    ## constraints ##
    constridx = 1

    # normalization constraint: E(1) = 1
    one = (Symbol=>Int64)[] :: SoSMonom
    #              idx    blk row  col  val
    setcon!(sdp,constridx, 1, one, one, 1.0)
    setrhs!(sdp,constridx, 1.0)
    constridx += 1

    # symmetry constraints
    for md in 0:d
        for monom in monoms(sos,md)
            dec = decomp(monom, div(d,2))
            (a1,b1) = dec[1]    # same as decomp1(monom,div(d,2))
            for (a,b) in dec
                if (a,b) == (a1,b1) || (a,b) == (b1,a1) # SDP solver already enforces matrix symmetry
                    continue
                end
                #@printf("symmetry: %s is %s\n", (a,b), (a1,b1))
                setcon!(sdp,constridx, 1, a, b, 1.0)
                setcon!(sdp,constridx, 1,a1,b1,-1.0)
                setrhs!(sdp,constridx,0.0)
                constridx += 1
            end
        end
    end

    # given constraints
    for poly in sos.constraints

        pd = deg(poly)
        if(pd > d)
            @printf("Warning: a constraint has degree %d, larger than solution degree %d",pd,d)
        end

        # promote to higher degree in all possible ways
        for promdeg in 0:(d-pd)
            for monom in monoms(sos,promdeg)

                for (k,v) in poly
                    (k1,k2) = decomp1(k * monom,div(d,2))
                    #                     blk row col val
                    setcon!(sdp, constridx, 1, k1, k2, v)
                end
                setrhs!(sdp, constridx, 0.0)
                constridx += 1
            end
        end
    end

    nvars = count_vars(sos, div(d,2))
    @printf("SDP with %d^2 entries and %d constraints\n", nvars, constridx-1)
    #show(sdp.cons)

    # solve
    sol = solve(sdp, solverinst)

    # TODO package this in some way
    println("Objective:")
    println(obj(sol))
    println("Primal:")
    println(primalmatrix(sol))
    println("Dual:")
    println(dualmatrix(sol))
end


function count_vars(sos :: SoS, deg :: Int64)
    n = length(sos.vars)
    sum([ binomial(d+n-1, d) for d in 0:deg])
end

