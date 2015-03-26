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
        mond = vcat( [ monoms(prog,i) for i in 0:d ] ... )
        mondrev = [ mond[i]=>i for i in 1:length(mond) ]
        dec = [ decomp(m, div(d,2)) for m in mond ]
    end


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
                for monom in monoms(prog,promdeg) # XXX should we avoid re-generating these for each constraint?
                    for (k,v) in poly
                        push!(CI, constridx)
                        push!(CJ, mondrev[k * monom])
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
        O = zeros(length(mond))
        for (k,v) in prog.objective
            O[mondrev[k]] = v
        end

        # separate the constant column from the rest
        Cconst = full(C[:,1])
        C = full(C[:,2:end])
        Oconst = O[1]
        O = O[2:end]
        mond = mond[2:end] # excluding constant term
        mondrev = [ mond[i]=>i for i in 1:length(mond) ]
        dec = dec[2:end]
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
        for j in 1:length(initial)
            for (a,b) in dec[j]
                setobj!(sdp, 1, a, b, initial[j])
            end
        end
    end


    # compute the nullspace
    @printf("Computing nullspace: ")
    @time B = null(C)
#    println(rref(C))


    # dual constraints = primal building-blocks
    @printf("Writing ideal basis -- ")
    @time for i in 1:size(B,2)
        
        # rewrite B[:,i] as a matrix in d/2 x d/2 moment decompositions
        for j in 1:size(B,1)

            for (a,b) in dec[j]
                setcon!(sdp, i, 1, a, b, B[j,i])
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
    
    CI = CJ = CV = C = Cconst = O = Oconst = B = initial = 0 # free some memory
    @time sol = solve(sdp, solverinst)


    # build & return a solution object
    @printf("Building solution object -- ")
    @time begin
        moments = [ mond[j] => sol.dualmatrix[1, dec[j][1]...] for j in 1:length(mond) ]
        moments[one] = 1.0
        SoSSolution(prog, d, sol.obj + adjust_obj, sol.dualobj + adjust_obj, moments, sol.primalmatrix)
    end
end
    


    #=
    @printf("computing rref: ")
    @time R = rref(C)
    Rnz = 0
    Rrows = 0
    for i in 1:size(R,2)
        v = R[:,i]

        ct = 0
        for val in v
            if abs(val) > 1e-7
                ct += 1
            end
        end

        if ct == 1
            continue
        end

        Rnz += ct + 1
        Rrows += 1
    end
    @printf("rref-nullspace has %d vectors, with %d nonzero entries (density %f)\n", Rrows, Rnz, Rnz / (Rrows * size(R,2)))
    =#


function count_monoms(prog :: Program, d :: Int64)
    n = length(prog.vars)
    binomial(d+n,d)
end

