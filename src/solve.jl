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


    @printf("Enumerating monomials -- ")
    @time begin
        # choose an order for indexing non-1 monomials
        monall = [ monoms(prog,i) for i in 0:d ]
        mon0d2 = vcat( monall[ 1:(div(d,2)+1) ]... ) # monomials of degree from 0 to d/2
    end

    
    @printf("Symmetrizing monomials -- ")
    @time omap, revomap = monom_orbits(vcat(monall...), prog.symmetries)


    @printf("Writing constraint-by-orbit matrix -- ")
    @time begin
        # write the (sparse) constraint-by-orbit occurence matrix C
        rows = Dict{Int64,Float64}[]
        for poly in prog.constraints
            pd = deg(poly)
            if(pd > d)
                @printf("Warning: a constraint has degree %d, larger than solution degree %d",pd,d)
            end

            # promote to higher degree in all possible ways
            for promdeg in 0:(d-pd)
                for monom in monall[promdeg+1] # +1 because monall is an array indexed from 1

                    row = Dict{Int64,Float64}()
                    for (k,v) in poly
                        key = omap[k*monom]
                        row[key] = v + get(row,key,0.0)
                    end

                    # only add this row if it's actually new
                    # XXX this is currently a major source of inefficiency! can we improve
                    # this without resorting to a full nearest-neighbor data structure?
                    # maybe hash constraints, exploiting the fact that duplicated rows will
                    # probably be bit-equal, not just epsilon-close, and that we don't have
                    # to be perfect...
                    if !any(o -> rowdist(row,o) < 1e-8, rows)
                        push!(rows,row)
                    end

                end
            end
        end

        # transcribe the rows to a matrix
        # XXX we have to specify width `length(revomap)` in case higher moments /
        # orbits don't show up in the promoted constraints. In principle,
        # though, it'd make the linear algebra faster to exclude these columns,
        # and manually re-include them after the linear algebra. At some point we
        # should do this; we can also reduce columns by joining orbits by
        # equality constraints (e.g. boolean constraints)
        # XXX we're not getting any mileage out of the sparsity here. Fix this.
        # We should be able to at least find `initial` using sparse QR, even if
        # we can't find the nullspace in a sparse-friendly way.
        C = zeros(length(rows),length(revomap))
        for i in 1:length(rows)
            for (j,v) in rows[i]
                C[i,j] = v
            end
        end
    end
    @printf("Size of C is %s\n", size(C))


    @printf("Manipulating data structures -- ")
    @time begin
        # write our objective in moments
        O = zeros(length(revomap))
        for (k,v) in prog.objective
            O[omap[k]] += v
        end

        # separate the constant column from the rest
#        Cconst = full(C[:,1])
#        C = full(C[:,2:end])
        Cconst = C[:,1]
        C = C[:,2:end]
        Oconst = O[1]
        O = O[2:end]
        # XXX omap
        revomap = revomap[2:end]
    end


    @printf("Precomputing decomposition-to-orbit map -- ")
    @time begin
        domap = Dict{(Int64,Int64),Int64}()
        for a in 1:length(mon0d2)
            for b in max(a,2):length(mon0d2)
                domap[(a,b)] = omap[mon0d2[a]*mon0d2[b]]-1
            end
        end
    end

    
    
    # find an initial solution for moments (not necessarily PSD)
    @printf("Computing initial value -- ")
    @time if size(C,1) > 0
        initial = vec(C \ Cconst)
    else
        initial = zeros(size(C,2))
    end


    # compute the nullspace
    @printf("Computing ideal basis -- ")
    @time B = null(C)


    # we will need to adjust our objective by the following constant term, which is missing from the SDP
    adjust_obj = Oconst - dot(O, initial)

    # set up problem & solver instance
    sdp = SDPSession(prog.maximize, size(B,2), length(mon0d2), adjust_obj)


    # write objective value
    @printf("Writing objective -- ")
    @time for i in 1:size(B,2)
        sdp_obj!(sdp, dot(O, B[:,i]))
    end


    # dual objective = primal constant-term
    # XXX TMP
    big = 0
    tot = 0

    @printf("Writing initial value -- ")
    @time begin
        sdp_con!(sdp, 0, 1, 1, -1.0)
        for ((a,b),orbit) in domap
            val = initial[orbit]
            tot += 1
            if(abs(val) > 1e-8)
                big += 1
                sdp_con!(sdp, 0, a, b, val)
            end
        end
    end
    @printf("Initial value: wrote %d of %d, density %f\n", big, tot, big/tot)


    # dual constraints = primal building-blocks
    # XXX TMP
    big = 0
    tot = 0

    @printf("Writing ideal basis -- ")
    @time for i in 1:size(B,2)
        
        # rewrite B[:,i] as a matrix in d/2 x d/2 moment decompositions
        for ((a,b),orbit) in domap
            val = B[orbit,i]
            tot += 1
            if(abs(val) > 1e-8)
                big += 1
                sdp_con!(sdp, i, a, b, val)
            end
        end

        @printf("done %d/%d constraints. wrote %d of %d, density %f\n", i, size(B,2), big, tot, big/tot)
    end


    # free some memory
    C = Cconst = O = Oconst = B = initial = monall = mond = mondrev = omap = revomap = domap = 0 # free some memory


    # solve the SDP
    @printf("Solving SDP with %d^2 entries and %d constraints... ", sdp.nmoments, sdp.nconstraints)
    @time sol = sdp_solve(sdp)


    # build & return a solution object
    @printf("Building solution object -- ")
    @time begin
        moments = Dict{SoSMonom,Float64}()
        for i in 1:length(mon0d2)
            for j in i:length(mon0d2)
                moments[mon0d2[i]*mon0d2[j]] = sol.primalmatrix[i,j]
            end
        end
        SoSSolution(prog, d, sol.primalobj, sol.dualobj, moments, sol.primalmatrix, sol.dualmatrix)
    end
end


rowdist(r1,r2) = halfrowdist(r1,r2) + halfrowdist(r2,r1)
function halfrowdist(r1 :: Dict{Int64,Float64}, r2 :: Dict{Int64,Float64})
    dist = 0.0
    for (k,v) in r1
        dist += (v - get(r2,k,0.0))^2
    end
    return dist
end

