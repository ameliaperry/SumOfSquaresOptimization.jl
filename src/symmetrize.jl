# symmetrization

function act(mon :: SoSMonom, perm :: Dict{Symbol,Symbol})
    [ get(perm,k,k) => v for (k,v) in mon ]
end

function act(poly :: SoSPoly, perm :: Dict{Symbol,Symbol})
    [ act(k,perm) => v for (k,v) in poly ]
end

function symmetrize!(sol :: SoSSolution, genperms :: Array{Dict{Symbol,Symbol}})

    # we need to do the following things:
    #  symmetrize the primal matrix -- over the diagonal action on pairs of degree <=d/2 monomials
    #  symmetrize the dual matrix -- over the same. (do these two together)
    #  symmetrize the dual vector "alpha" entries -- actually this is impossible to do, since the set of
    #    alpha constraints / matrices isn't symmetric, and can't be chosen symmetrically without heavy
    #    redundancy. (think of degree 4 SoS, the monomial x1 x3 x5 under cyclic symmetry of x1 ... x6 --
    #    we can't choose a single decomposition (decomp1) of this into a pair of degree 2 monomials in a
    #    symmetric way.)
    #    For now, we leave the alpha vectors in an inconsistent state! XXX TODO, better hope nobody cares.
    #    What we should probably do is at least make them consistent again, even if they can't be made
    #    symmetric. Does this require solving a linear system...? Is there a clean way?
    #  symmetrize the dual vector -- have to act on constraint polynomials and their promoting monomials.
    
    
    # symmetrize primal & dual matrix
    # iterate over pairs of monomials. if we haven't seen a pair, symmetrize its orbit

    halfmonoms = [ SumOfSquaresOptimization.monoms(sol.prog,i) for i in 0:div(sol.degree,2) ]
    halfmonoms = vcat(halfmonoms...)
    
    seen = Set{(SoSMonom,SoSMonom)}()
    for b1 in halfmonoms
        for b2 in halfmonoms
            if in((b1,b2), seen)
                continue
            end
            
            orbit = Set{(SoSMonom,SoSMonom)}()
            q = (SoSMonom,SoSMonom)[]
            primaltot = 0.0
            dualtot = 0.0
            
            # start off with b1 and b2
            push!(orbit,(b1,b2))
            push!(q,(b1,b2))
            
            # explore the orbit
            while !isempty(q)
                (c1,c2) = shift!(q)
                
                # note the primal & dual values here
                primaltot += primalmatrix(sol.sdp)[1][c1,c2]
                dualtot += dualmatrix(sol.sdp)[1][c1,c2]
                
                # act every possible way on the current element
                for perm in genperms
                    (i1,i2) = (act(c1,perm), act(c2,perm))
                    
                    if !in((i1,i2),orbit)
                        push!(orbit,(i1,i2))
                        push!(q,(i1,i2))
                    end
                end
            end
            
            # now set the matrix values to the average
            primalavg = primaltot / length(orbit)
            dualavg = dualtot / length(orbit)
            
            for (c1,c2) in orbit
                primalmatrix(sol.sdp)[1][c1,c2] = primalavg
                dualmatrix(sol.sdp)[1][c1,c2] = dualavg
            end
            
            union!(seen,orbit)
        end
    end
    

    # symmetrize the "alpha" entries of the dual vector
    #
    # XXX this doesn't work as written below: when we act on
    # (a,b,c,d) by the diagonal action, we don't necessarily get
    # another quadruple for which an alpha variable exists. If we
    # included alpha variables for every quadruple, the SDP would
    # become huge. We can't choose a minimal subset in a symmetric
    # way.
    #=
    seen = Set{(SoSMonom,SoSMonom,SoSMonom,SoSMonom)}()
    for ((ba,bb,bc,bd),bi) in sol.alphaidx
        if in((ba,bb,bc,bd), seen)
            continue
        end

        orbit = Set{(SoSMonom,SoSMonom,SoSMonom,SoSMonom)}()
        q = (SoSMonom,SoSMonom,SoSMonom,SoSMonom)[]
        tot = 0.0
        
        # start off with b1 and b2
        push!(orbit,(ba,bb,bc,bd))
        push!(q,(ba,bb,bc,bd))
        
        # explore the orbit
        while !isempty(q)
            (ca,cb,cc,cd) = shift!(q)
            
            # note the alpha value here
            tot += dualvector(sol.sdp)[sol.alphaidx[(ca,cb,cc,cd)]]
            
            # act every possible way on the current element
            for perm in genperms
                (ia,ib,ic,id) = (act(ca,perm), act(cb,perm), act(cc,perm), act(cd,perm))
                
                if !in((ia,ib,ic,id),orbit)
                    push!(orbit,(ia,ib,ic,id))
                    push!(q,(ia,ib,ic,id))
                end
            end
        end
        
        # now set the matrix values to the average
        avg = tot / length(orbit)
        
        for (ca,cb,cc,cd) in orbit
            dualvector(sol.sdp)[sol.alphaidx[(ca,cb,cc,cd)]] = avg
        end
        
        union!(seen,orbit)
    end
    =#
     


    # symmetrize the "beta" entries of the dual vector
    seen = Set{(SoSPoly,SoSMonom)}()
    for ((bp,bm),bi) in sol.betaidx
        if in((bp,bm), seen)
            continue
        end

        orbit = Set{(SoSPoly,SoSMonom)}()
        q = (SoSPoly,SoSMonom)[]
        tot = 0.0
        
        # start off with b1 and b2
        push!(orbit,(bp,bm))
        push!(q,(bp,bm))
        
        # explore the orbit
        while !isempty(q)
            (cp,cm) = shift!(q)
            
            # note the primal & dual values here
            if !haskey(sol.betaidx, (cp,cm))
                @printf("Error symmetrizing dual vector.")
                @printf("Symmetry violation: %s is a constraint polynomial, but not %s\n", bp, cp)
                @printf("SoS solution object is now in an inconsistent state.")
                throw(ArgumentError("SoS program does not have these symmetries"))
            end
            tot += dualvector(sol.sdp)[sol.betaidx[(cp,cm)]]
            
            # act every possible way on the current element
            for perm in genperms
                (ip,imo) = (act(cp,perm), act(cm,perm))
                
                if !in((ip,imo),orbit)
                    push!(orbit,(ip,imo))
                    push!(q,(ip,imo))
                end
            end
        end
        
        # now set the matrix values to the average
        avg = tot / length(orbit)
        
        for (cp,cm) in orbit
            dualvector(sol.sdp)[sol.betaidx[(cp,cm)]] = avg
        end
        
        union!(seen,orbit)
    end

end

function symmetrize_dihedral!(sol :: SoSSolution, cycle :: Array{Symbol})
    n = length(cycle)
    rotate = (Symbol=>Symbol)[ cycle[i] => cycle[(i % n) + 1] for i in 1:n ]
    flip = (Symbol=>Symbol)[ cycle[i] => cycle[n + 1 - i] for i in 1:n ]
    symmetrize!(sol, [rotate,flip])
end

function symmetrize_cyclic!(sol :: SoSSolution, cycle :: Array{Symbol})
    n = length(cycle)
    rotate = (Symbol=>Symbol)[ cycle[i] => cycle[(i % n) + 1] for i in 1:n ]
    symmetrize!(sol, [rotate])
end

function symmetrize_full!(sol :: SoSSolution, syms :: Array{Symbol})
    if(length(syms) < 2)
        return
    end

    n = length(cycle)
    rotate = (Symbol=>Symbol)[ syms[i] => syms[(i % n) + 1] for i in 1:n ]
    transpose = (Symbol=>Symbol){ syms[1] => syms[2], syms[2] => syms[1] }
    symmetrize!(sol, [rotate,transpose])
end

