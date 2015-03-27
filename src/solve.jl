#####################
###### Solving ######
#####################


#  Make the magic happen.
function sossolve(prog :: Program, d :: Int64; solver="csdp")

    # sanity check
    if(d <= 0)
        throw(ArgumentError("Degree must be positive"))
    elseif(d < deg(prog.objective))
        throw(ArgumentError("Degree must be at least that of the objective"))
    elseif(d % 2 == 1)
        throw(ArgumentError("Degree must be even"))
    elseif(length(prog.vars) == 0)
        throw(ArgumentError("Program is empty"))
    end


    @printf("Enumerating monomials and their decompositions -- ")
    @time begin
        # choose an order for indexing non-1 monomials
        monall = [ monoms(prog,i) for i in 0:d ]
        mond = vcat(monall...)
        mondrev = [ mond[i]=>i for i in 1:length(mond) ]
        mon0d2 = vcat( monall[ 1:(div(d,2)+1) ]... ) # monomials of degree from 1 to d/2
        mon1d2 = mon0d2[2:end]
        dec = [ decomp(m, div(d,2)) for m in mond ] # XXX this is only used in one place now. can we push it out?
    end


    @printf("Symmetrizing monomials -- ")
    @time omap, revomap = monom_orbits(mond, prog.symmetries)


    @printf("Writing constraint-by-moment matrix -- ")
    @time begin
        # write the (sparse) constraint-by-moment occurence matrix C
        CI = Int64[]
        CJ = Int64[]
        CV = Float64[]
        constridx=1
        for poly in prog.constraints
            pd = deg(poly)
            if(pd > d)
                @printf("Warning: a constraint has degree %d, larger than solution degree %d",pd,d)
            end

            # promote to higher degree in all possible ways
            for promdeg in 0:(d-pd)
                for monom in monall[promdeg+1] # +1 because monall is an array indexed from 1
                    for (k,v) in poly
                        push!(CI, constridx)
                        push!(CJ, omap[k * monom])
                        push!(CV, v)
                    end
                    constridx += 1
                end
            end
        end
        C = sparse(CI,CJ,CV)
    end


    @printf("Manipulating data structures -- ")
    @time begin
        # write our objective in moments
        O = zeros(length(revomap))
        for (k,v) in prog.objective
            O[omap[k]] += v
        end

        # separate the constant column from the rest
        Cconst = full(C[:,1])
        C = full(C[:,2:end])
        Oconst = O[1]
        O = O[2:end]
        mond = mond[2:end] # excluding constant term
        mondrev = [ mond[i]=>i for i in 1:length(mond) ] # XXX unused
        dec = dec[2:end]
        # XXX omap
        revomap = revomap[2:end]
    end


    # set up problem & solver instance
    sdp = SparseSDP(maximize = !prog.maximize) # negate maximization, because our primal is the SDP dual
    if (solver == "sdpa")
        solverinst = SDPA()
    elseif (solver == "csdp")
        solverinst = CSDP()
    else
        throw(ArgumentError("Unsupported solver."))
    end
    
    
    # find an initial solution for moments (not necessarily PSD)
    @printf("Computing initial value: ")
    @time initial = vec(C \ Cconst)

    # dual objective = primal constant-term
    @printf("Writing initial value: ")
    @time begin
        setobj!(sdp, 1, one, one, -1.0)
        for a in mon0d2
            for b in mon1d2
                setobj!(sdp, 1, a, b, initial[omap[a*b]-1])
            end
        end
    end


    # compute the nullspace
    @printf("Computing ideal basis: ")
    @time B = null(C)
#    println(rref(C))


    # dual constraints = primal building-blocks
    @printf("Writing ideal basis -- ")
    @time for i in 1:size(B,2)
        
        # rewrite B[:,i] as a matrix in d/2 x d/2 moment decompositions
        for a in mon0d2
            for b in mon1d2
                setcon!(sdp, i, 1, a, b, B[omap[a*b]-1, i])
            end
        end


        rhs = dot(O, B[:,i])
        setrhs!(sdp, i, rhs)
    end


    # we will need to adjust our objective by the following constant term, which is missing from the SDP
    adjust_obj = Oconst - dot(O, initial)


    # solve the SDP
    nvars = count_monoms(prog, div(d,2))
    @printf("Solving SDP with %d^2 entries and %d constraints... ", count_monoms(prog, div(d,2)), size(B,2))
    
    CI = CJ = CV = C = Cconst = O = Oconst = B = initial = monall = mondrev = omap = revomap = 0 # free some memory
    @time sol = solve(sdp, solverinst)


    # build & return a solution object
    @printf("Building solution object -- ")
    @time begin
        moments = [ mond[j] => sol.dualmatrix[1, dec[j][1]...] for j in 1:length(mond) ]
        moments[one] = 1.0
        SoSSolution(prog, d, sol.obj + adjust_obj, sol.dualobj + adjust_obj, moments, sol.primalmatrix)
    end
end


function count_monoms(prog :: Program, d :: Int64)
    n = length(prog.vars)
    binomial(d+n,d)
end

