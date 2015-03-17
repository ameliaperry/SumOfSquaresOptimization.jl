#######################################
###### Monomials and polynomials ######
#######################################

#  An SoSMonom represents a monomial: it maps each symbol to its exponent.
#  Exponents should always be strictly positive (never 0).
typealias SoSMonom Dict{Symbol,Int64}

#  An SoSPoly maps each monomial to its coefficient.
typealias SoSPoly Dict{SoSMonom,Float64}


#  Multiplication of polynomials.
#    there should be a faster sum here that manages iterators, a la mergesort
function *(a::SoSMonom, b::SoSMonom)
    # A, including overlap with B
    map = [k => (v + get(b,k,0)) for (k,v) in a] :: SoSMonom
    # B \ A
    for (k,v) in b
        if(! haskey(a,k))
            map[k] = v
        end
    end
    return map
end

#  Addition of polynomials
#    there should be a faster sum here that manages iterators, a la mergesort
#    also, this is the same code as multiplication of SosMonoms -- make non-redundant
function +(a::SoSPoly, b::SoSPoly)
    # A, including overlap with B
    map = [k => (v + get(b,k,0)) for (k,v) in a] :: SoSPoly
    # B \ A
    for (k,v) in b
        if(! haskey(a,k))
            map[k] = v
        end
    end
    return map
end

#  Addition in-place of polynomials
function addpoly!(a::SoSPoly, b::SoSPoly)
    for (k,v) in b
        a[k] = get(a,k,0.0) + v
    end
end


#  Multiplication of polynomials
#    just the naive method
function *(a::SoSPoly, b::SoSPoly)
    # there are maybe more efficient ways to do this
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
function deg(m :: SoSMonom)
    sum(values(m))
end
function deg(p :: SoSPoly)
    maximum([deg(k) for (k,v) in p])
end


#  Enumerate all monomials of degree d.
#    Julia knows about combinations (without replacement), and we
#    want to choose variables with replacement. We simulate this
#    via (v+d-1 choose d).
function monoms(varr :: Array{Symbol}, d :: Int64)
    v = length(varr)
    [combToMonom(varr,co) for co in combinations([1:(v+d-1)],d)]
end
function monoms(d :: Int64) # edge case
    if (d == 0)
        return [(Symbol=>Int64)[]]
    else
        return []
    end
end
function combToMonom(varr :: Array{Symbol}, co :: Array{Int64})
    m = (Symbol=>Int64)[] :: SoSMonom
    shift = 0
    for i in co
        k = varr[i+shift]
        m[k] = 1 + get(m,k,0);
        shift -= 1;
    end
    m
end



#  Find all decompositions of mon into two monomials (m1,m2)
#    each of degree at most degmax.
#    TODO: implement as an iterator.
function decomp(mon :: SoSMonom, degmax :: Int64)
    # maybe there's room for improvement here ...
    list = (SoSMonom,SoSMonom)[]
    varr = [k for (k,v) in mon]
    
    # for each candidate m1, see if its complement m2 is sensible
    mondeg = deg(mon)
    dmin = max(0,mondeg-degmax)
    dmax = min(mondeg,degmax)
    for deg in dmin:dmax
        for m1 in monoms(varr,deg)

            # construct the complement
            m2 = [k=>(v-get(m1,k,0)) for (k,v) in mon]

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

#  Find a single decomposition of mon into two monomials (m1,m2)
#    each of degree at most degmax.
#    TODO: make this code less redundant. Maybe once decomp() is
#    an iterator, we can just pull off the first one.
function decomp1(mon :: SoSMonom, degmax :: Int64)
    # maybe there's room for improvement here ...
    varr = [k for (k,v) in mon]
    
    # for each candidate m1, see if its complement m2 is sensible
    mondeg = deg(mon)
    dmin = max(0,mondeg-degmax)
    dmax = min(mondeg,degmax)
    for deg in dmin:dmax
        for m1 in monoms(varr,deg)

            # construct the complement
            m2 = [k=>(v-get(m1,k,0)) for (k,v) in mon]

            # test for being sensible (ie non-negative exponents), and remove any zero entries
            if(length(m2) > 0 && minimum(values(m2)) < 0)
                continue
            end

            # remove any zero entries
            filter!((k,v)->(v!=0), m2)
            
            return (m1,m2)
        end
    end

    throw(Error("No monomial decompositions found! Some degree must be off. Please report this bug."))
end

