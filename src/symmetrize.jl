# symmetrization

function act(mon :: SoSMonom, perm :: Dict{Symbol,Symbol})
    [ get(perm,k,k) => v for (k,v) in mon ] :: SoSMonom
end

# XXX currently unused
function act(poly :: SoSPoly, perm :: Dict{Symbol,Symbol})
    [ act(k,perm) => v for (k,v) in poly ] :: SoSPoly
end

# returns a Dict{SoSMonom,Symbol} mapping monomials to their orbits
function monom_orbits(prog :: Program, genperms :: Array{Dict{Symbol,Symbol}}, deg :: Int64)

    monoms = vcat( [ SumOfSquaresOptimization.monoms(prog,i) for i in 0:deg ] ... )
    ret = Dict{SoSMonom,Symbol}()
    
    seen = Set{SoSMonom}()
    for b in monoms
        if in(b, seen)
            continue
        end
        
        orbit = Set{SoSMonom}()
        q = SoSMonom[]
        
        # start off with b
        push!(orbit, b)
        push!(q, b)
        
        # explore the orbit
        while !isempty(q)
            c = shift!(q)
            
            # act every possible way on the current element
            for perm in genperms
                i = act(c,perm)
                
                if !in(i,orbit)
                    push!(orbit,i)
                    push!(q,i)
                end
            end
        end

        sym = gensym()
        for i in orbit
            ret[i] = sym
        end
        
        union!(seen,orbit)
    end
    
    ret
end


function symmetrize!(prog :: Program, genperms :: Array{Dict{Symbol,Symbol}})
    # TODO
end


function symmetrize_dihedral!(prog :: Program, cycle :: Array{Symbol})
    n = length(cycle)
    rotate = (Symbol=>Symbol)[ cycle[i] => cycle[(i % n) + 1] for i in 1:n ]
    flip = (Symbol=>Symbol)[ cycle[i] => cycle[n + 1 - i] for i in 1:n ]
    symmetrize!(prog, [rotate,flip])
end

function symmetrize_cyclic!(prog :: Program, cycle :: Array{Symbol})
    n = length(cycle)
    rotate = (Symbol=>Symbol)[ cycle[i] => cycle[(i % n) + 1] for i in 1:n ]
    symmetrize!(prog, [rotate])
end

function symmetrize_full!(prog :: Program, syms :: Array{Symbol})
    if(length(syms) < 2)
        return
    end

    n = length(cycle)
    rotate = (Symbol=>Symbol)[ syms[i] => syms[(i % n) + 1] for i in 1:n ]
    transpose = (Symbol=>Symbol){ syms[1] => syms[2], syms[2] => syms[1] }
    symmetrize!(prog, [rotate,transpose])
end

