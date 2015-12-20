#####################
###### Solving ######
#####################

sosdebug = false

macro timeifdebug(str,ex)
    if sosdebug
        quote
            @printf("%s -- ", $(esc(str)))
            @time $(esc(ex))
        end
    else
        quote $(esc(ex)) end
    end
end

function row_redundancy_object(wid)
    randproj = rand(1,wid)
    rowsd = SortedDict(Dict{Float64,SparseMatrixCSC{Float64,Int64}}())
    (randproj, rowsd)
end

function row_redundancy_check(redobj, row)
    (randproj, rowsd) = redobj
    rval = (randproj * row)[1,1]
    
    st = searchsortedfirst(rowsd,rval)
    while status((rowsd,st)) == 1
        (k,v) = deref((rowsd,st))
        if k - rval > 1e-7
            break
        end
        if norm(v - row) < 1e-7
            return false
        end
        st = advance((rowsd,st))
    end
    
    st = searchsortedlast(rowsd,rval)
    while status((rowsd,st)) == 1
        (k,v) = deref((rowsd,st))
        if rval - k > 1e-7
            break
        end
        if norm(v - row) < 1e-7
            return false
        end
        st = regress((rowsd,st))
    end

    rowsd[rval] = row
    true
end


#  Make the magic happen.
function sossolve(prog :: Program, d :: Int64; solver="csdp", call="csdp")

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


    @timeifdebug "Enumerating monomials" begin
        # choose an order for indexing non-1 monomials
        monall = [ monoms(prog,i) for i in 0:d ]
        mon0d2 = vcat( monall[ 1:(div(d,2)+1) ]... ) # monomials of degree from 0 to d/2
    end

    
    @timeifdebug "Symmetrizing monomials" begin
        omap, revomap = monom_orbits(vcat(monall...), prog.symmetries)
    end


    @timeifdebug "Writing constraint-by-orbit matrix" begin
        # choose a random projection for testing row redundancy

        # write the (sparse) constraint-by-orbit occurence matrix C, removing redundant constraints
        CI = Array{Int64,1}()
        CJ = Array{Int64,1}()
        CV = Array{Float64,1}()
        redund = row_redundancy_object(length(revomap))
        i = 1
        for poly in prog.constraints
            pd = deg(poly)
            if pd > d
                println("Warning: a constraint has degree $pd, larger than relaxation degree $d")
            end

            # promote to higher degree in all possible ways
            for promdeg in 0:(d-pd)
                for monom in monall[promdeg+1] # +1 because monall is an array indexed from 1

                    row = spzeros(length(revomap),1)
                    for (k,v) in poly
                        key = omap[k*monom]
                        row[key,1] = v + get(row,key,0.0)
                    end

                    # only add this row if it's actually new
                    if row_redundancy_check(redund, row)
                        rv = rowvals(row)
                        append!(CI,fill(i,length(rv)))
                        append!(CJ,rv)
                        append!(CV,nonzeros(row))
                        i += 1
                    end

                end
            end
        end

        # XXX we have to specify width `length(revomap)` in case higher moments /
        # orbits don't show up in the promoted constraints. In principle,
        # though, it'd make the linear algebra faster to exclude these columns,
        # and manually re-include them after the linear algebra. At some point we
        # should do this; we can also reduce columns by joining orbits by
        # equality constraints (e.g. boolean constraints)
        C = sparse(CI,CJ,CV,i-1,length(revomap))
    end
    @printf("Size of C is %s\n", size(C))


    @timeifdebug "Manipulating data structures" begin
        # write our objective in moments
        O = zeros(length(revomap))
        for (k,v) in prog.objective
            O[omap[k]] += v
        end

        # separate the constant column from the rest
        Cconst = C[:,1]
        C = C[:,2:end]
        Oconst = O[1]
        O = O[2:end]
        # XXX omap
        revomap = revomap[2:end]
    end


    @timeifdebug "Precomputing decomposition-to-orbit map" begin
        domap = Dict{Tuple{Int64,Int64},Int64}()
        for a in 1:length(mon0d2)
            for b in max(a,2):length(mon0d2)
                domap[(a,b)] = omap[mon0d2[a]*mon0d2[b]]-1
            end
        end
    end

    
    
    # find an initial solution for moments (not necessarily PSD)
    @timeifdebug "Computing initial value" if size(C,1) > 0
        initial = vec(C \ full(Cconst))
        # test that this is actually a solution
        if norm(Cconst - C * initial) > 1e-6
            #  XXX TODO
        end
    else
        initial = zeros(size(C,2))
    end


    # compute the nullspace
    @timeifdebug "Computing ideal basis" begin
        B = nullspace(full(C))
    end


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
    println("Initial value: density $(big/tot)")


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

        println("done $i/$(size(B,2)) constraints, density $(big/tot)")
    end


    # free some memory
    C = Cconst = O = Oconst = B = initial = monall = mond = mondrev = omap = revomap = domap = 0


    # solve the SDP
    @timeifdebug "Solving SDP with $(sdp.nmoments)^2 entries and $(sdp.nconstraints) constraints" begin
        sol = sdp_solve(sdp,call=call)
    end


    # build & return a solution object
    @timeifdebug "Building solution object" begin
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

