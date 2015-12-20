#######################################
###### Monomials and polynomials ######
#######################################

#  An SoSMonom represents a monomial: it maps each symbol to its exponent.
#  Exponents should always be strictly positive (never 0).
typealias SoSMonom Dict{Symbol,Int64}

#  An SoSPoly maps each monomial to its coefficient.
typealias SoSPoly Dict{SoSMonom,Float64}

#    there should be a faster option here that manages iterators, a la mergesort
macro dictplus(a, b)
    quote
        # A, including overlap with B
        map = Dict([ k => v + get(b,k,0) for (k,v) in a ])
        # B \ A
        for (k,v) in b
            if !haskey(a,k)
                map[k] = v
            end
        end
        map
    end
end

#  Multiplication of monomials.
*(a::SoSMonom,b::SoSMonom) = @dictplus(a,b)
+(a::SoSPoly,b::SoSPoly) = @dictplus(a,b)

#  Addition in-place of polynomials
function addpoly!(a::SoSPoly, b::SoSPoly)
    for (k,v) in b
        a[k] = get(a,k,0.0) + v
    end
end

#  Multiplication of polynomials -- just the naive method
function *(a::SoSPoly, b::SoSPoly)
    map = SoSPoly()
    for (k1,v1) in a
        for (k2,v2) in b
            kp = k1 * k2
            map[kp] = get(map,kp,0) + v1 * v2
        end
    end
    return map
end

#  Degree of monomial / polynomial
deg(m :: SoSMonom) = sum(values(m))
deg(p :: SoSPoly) = (length(p) == 0) ? 0 : maximum([deg(k) for (k,v) in p]) 


#  Enumerate all monomials of degree d.
#    Julia knows about combinations (without replacement), and we
#    want to choose variables with replacement. We simulate this
#    via (v+d-1 choose d).
monoms(varr :: Array{Symbol}, d :: Int64) =
    [ combToMonom(varr,co) for co in combinations(1:(length(varr)+d-1),d) ] :: Array{SoSMonom,1}
monoms(d :: Int64) = (d == 0) ? [SoSMonom()] : Array{SoSMonom,1}() # edge case

function combToMonom(varr :: Array{Symbol}, co :: Array{Int64})
    m = Dict{Symbol,Int64}() :: SoSMonom
    shift = 0
    for i in co
        k = varr[i+shift]
        m[k] = 1 + get(m,k,0);
        shift -= 1;
    end
    m
end


# Find all decompositions of mon into two monomials (m1,m2), each of degree
# at most degmax.
# XXX currently unused
function decomp(mon :: SoSMonom, degmax :: Int64)
    # maybe there's room for improvement here ...
    list = Tuple{SoSMonom,SoSMonom}[]
    varr = [k for (k,v) in mon]
    
    # for each candidate m1, see if its complement m2 is sensible
    mondeg = deg(mon)
    dmin = max(0,mondeg-degmax)
    dmax = min(mondeg,degmax)
    for deg in dmin:dmax
        for m1 in monoms(varr,deg)

            # construct the complement
            m2 = Dict([ k => (v-get(m1,k,0)) for (k,v) in mon ])

            # test for being sensible (ie non-negative exponents), and remove any zero entries
            if(length(m2) > 0 && minimum(values(m2)) < 0)
                continue
            end

            # remove any zero entries
            filter!((k,v)->(v!=0), m2)
            
            push!(list, (m1,m2))
        end
    end

    list
end

