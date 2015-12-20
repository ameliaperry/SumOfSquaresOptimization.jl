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


# The row-redundancy test works as follows: we choose a random 1-d projection,
# and keep a dictionary of existing rows sorted by their value in that 1-d
# projection. Each time we get a new row, we take its value in that projection
# and fully compare that row with the projection-nearby rows. This should result in
# very few superfluous full-row comparisons.
# In order to get the nearby rows as indexed by projection value, we use the SortedDict
# type from DataStructures.jl.
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
                        orbit = omap[k*monom]
                        row[orbit,1] += v 
                    end

                    # add this row, but only if it's actually new
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

        # we have to specify width `length(revomap)` in case higher moments /
        # orbits don't show up in the promoted constraints. 
        # XXX In principle, though, it'd make the linear algebra faster to
        # exclude these columns, and manually re-include them after the linear
        # algebra. At some point we should do this; we can also reduce columns
        # by joining orbits by equality constraints (e.g. boolean constraints)
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
        revomap = revomap[2:end]
        
        # decomposition-to-orbit map
        domap = Dict{Tuple{Int64,Int64},Int64}()
        for a in 1:length(mon0d2)
            for b in max(a,2):length(mon0d2)
                domap[(a,b)] = omap[mon0d2[a]*mon0d2[b]]-1
            end
        end

        omap = 0 # done with this
    end
    
    
    # find an initial solution for moments -- satisfying affine constraints, but not necessarily PSD
    @timeifdebug "Computing initial value" if size(C,1) > 0
        initial = vec(C \ full(Cconst))
        # test that this actually satisfies constraints. if not, declare infeasibility
        if norm(Cconst - C * initial) > 1e-6

            #  XXX TODO
        end
    else
        initial = zeros(size(C,2))
    end


    # compute the nullspace
    # XXX really want to do this in a sparse way, from the QR decomposition of B
    @timeifdebug "Computing ideal basis" begin
        B = nullspace(full(C))
    end


    # we will need to adjust our objective by the following constant term, which is missing from the SDP
    adjust_obj = Oconst - dot(O, initial)

    # set up problem & solver instance
    sdp = SDPSession(prog.maximize, size(B,2), length(mon0d2), adjust_obj)


    # write objective value
    @timeifdebug "Writing objective -- " begin
        for i in 1:size(B,2)
            sdp_obj!(sdp, dot(O, B[:,i]))
        end
    end


    # dual objective = primal constant-term
    @timeifdebug "Writing initial value -- " begin
        sdp_con!(sdp, 0, 1, 1, -1.0)
        for ((a,b),orbit) in domap
            val = initial[orbit]
            if(abs(val) > 1e-8)
                sdp_con!(sdp, 0, a, b, val)
            end
        end
    end

    # dual constraints = primal building-blocks
    @timeifdebug "Writing ideal basis -- " begin
        for i in 1:size(B,2)
            # rewrite B[:,i] as a matrix in d/2 x d/2 moment decompositions
            for ((a,b),orbit) in domap
                val = B[orbit,i]
                if(abs(val) > 1e-8)
                    sdp_con!(sdp, i, a, b, val)
                end
            end
        end
    end


    # free some memory
    C = Cconst = O = Oconst = B = initial = monall = mond = mondrev = revomap = domap = 0


    # solve the SDP
    println("Solving SDP with $(sdp.nmoments)^2 entries and $(sdp.nconstraints) constraints") 
    @time sdpsol = sdp_solve(sdp,call=call)


    # build & return a solution object
    @timeifdebug "Building solution object" begin
        if sdpsol.status == :Infeasible
            return sos_sol_unbounded(prog,d)
        elseif sdpsol.status == :Unbounded
            return sos_sol_infeasible(prog,d)
        end
        moments = Dict{SoSMonom,Float64}()
        for i in 1:length(mon0d2)
            for j in i:length(mon0d2)
                moments[mon0d2[i]*mon0d2[j]] = sdpsol.primalmatrix[i,j]
            end
        end
        status = sdpsol.status == :Infeasible ? :Unbounded : sdpsol.status == :Unbounded ? :Infeasible : sdpsol.status
        SoSSolution(prog, d, sdpsol.primalobj, sdpsol.dualobj, moments, sdpsol.primalmatrix, sdpsol.dualmatrix, status)
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

