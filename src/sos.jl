#############################
###### Sums of Squares ######
#############################

#  A system of polynomial equations.
#    We convert inequalities into equations with slack variables (squared).
type SoS
    slackvars :: Int  # counter
    vars :: Set{Symbol}  # all variables, including slack variables
    constraints :: Array{SoSPoly}  # polynomials that are constrained to 0
    namedconstraints :: Dict{Any,SoSPoly}
    objective :: SoSPoly
    maximize :: Bool
end
SoS(; maximize=true, minimize=false) = begin
    # default to maximize; switch to minimize if there's any sign to do so
    mx = true
    if(maximize == false || minimize == true)
        mx = false
    end
    SoS(0, Set{Symbol}(), [], Dict{Any,SoSPoly}(), SoSPoly(), mx)
end


#  Generate a new slack variable (a Symbol)
function slackvar(sos :: SoS)
    sos.slackvars += 1
    # there's a better symbol() in julia 0.4, but let's target 0.3
    return symbol(@sprintf("slack%d",sos.slackvars))
end

#  Convenience: enumerate monomials in the variables that an SoS
#    program knows about.
function monoms(sos :: SoS, d :: Int64)
    varr = [v for v in sos.vars]
    monoms(varr,d)
end


typealias Maybe{T} Union(T,Nothing)




#  Parse a polynomial (in Expr form) into a SoSPoly, by traversing
#    a Julia syntax tree.
function parsepoly(sos :: Maybe{SoS}, ex::Symbol)
    if(sos != nothing)
        push!(sos.vars, ex)
    end

    return [[ex => 1] => 1.0] :: SoSPoly
end
function parsepoly{T<:Number}(sos :: Maybe{SoS}, ex :: T)
    return [(Symbol=>Int64)[] => convert(Float64, ex)] :: SoSPoly
end
function parsepoly(sos :: Maybe{SoS}, ex :: Expr)
    if ex.head == :call
        if(ex.args[1] == :(+))
            ret = parsepoly(sos,ex.args[2])
            for i in 3:length(ex.args)
                ret += parsepoly(sos,ex.args[i])
            end
            return ret
        elseif(ex.args[1] == :(*))
            ret = parsepoly(sos,ex.args[2])
            for i in 3:length(ex.args)
                ret *= parsepoly(sos,ex.args[i])
            end
            return ret
        elseif(ex.args[1] == :(-))
            if(length(ex.args) == 2) # unary minus
                return parsepoly(sos, :( (-1) * $(ex.args[2]) ))
            elseif(length(ex.args) == 3) # binary minus
                return parsepoly(sos, :( $(ex.args[2]) + (-1) * $(ex.args[3]) ))
            end
        elseif(ex.args[1] == :(^))
            if(! isa(ex.args[3],Int64) || ex.args[3] < 0)
                throw(ArgumentError(@sprintf("Bad exponent: %s. Exponents must be non-negative integers", args[3])))
            elseif(ex.args[3] == 0)
                return parsepoly(sos, 1)
            end
            operand = parsepoly(sos, ex.args[2])
            return prod([operand for i in 1:ex.args[3]])
        end
    elseif ex.head == :block # this happens on the RHS of single-equals
        if length(ex.args) > 2
            throw(ArgumentError(@sprintf("Multi-line block in expression: %s", ex)))
        end
        return parsepoly(sos, ex.args[2]) # assume it's a 1-line block
    end
    throw(ArgumentError(@sprintf("Expression not understood as polynomial: %s", ex)))
end


#  Massage any constraint (given in Expr form) into a polynomial equaling zero.
function relationToPolynomial(sos :: SoS, ex :: Expr)
    if ex.head == :comparison
        if(ex.args[2] == :(>=) || ex.args[2] == :(≥))
            polyex = :( $(ex.args[1]) - $(ex.args[3]) - $(slackvar(sos))^2 )
        elseif(ex.args[2] == :(==))
            polyex = :( $(ex.args[1]) - $(ex.args[3]) )
        elseif(ex.args[2] == :(<=) || ex.args[2] == :(≤))
            polyex = :( $(ex.args[3]) - $(ex.args[1]) - $(slackvar(sos))^2 )
        else
            throw(ArgumentError(
                "Comparison operator must be >=, ≥, ==, =, <=, or ≤"))
        end
    elseif ex.head == :(=) # isn't Julia comparison, but we should support
        polyex = :( $(ex.args[1]) - $(ex.args[2]) )
    else
        throw(ArgumentError("Constraint must be an equation or inequality"))
    end

    parsepoly(sos, polyex)
end



function constraint(sos,str)
    poly = relationToPolynomial(sos, parse(str))
    push!(sos.constraints, poly)
end

function partconstraint(sos,key,str)
    poly = parsepoly(sos, parse(str))
#    poly = relationToPolynomial(sos, parse(str))

    if haskey(sos.namedconstraints, key)
        prev = sos.namedconstraints[key]
        addpoly!(prev,poly)
    else
        push!(sos.constraints, poly)
        sos.namedconstraints[key] = poly
    end
end

function partobjective(sos,str)
    poly = parsepoly(sos, parse(str))
    sos.objective += poly
end

function objective(sos,str)
    poly = parsepoly(sos, parse(str))
    sos.objective = poly
end


function perturb_objective(sos,tol,deg)
    for md in 1:deg
        for mon in monoms(sos,md)
            poly = [ mon => tol * randn() ] :: SoSPoly
            sos.objective += poly
        end
    end
end

function perturb_objective_sparse(sos,tol,s,deg)
    mon = [ monoms(sos,i) for i in 0:deg ]
    mon = vcat(mon...)
    for i in 1:s
        idx = rand(1:length(mon))
        m = mon[idx]
        poly = [ m => tol * randn() ] :: SoSPoly
        sos.objective += poly
    end
end


#  The following four macros are the main user-facing ways to prepare
#    a SoS program.
#    This section is `esc` hell and was a pain to figure out.

#  Add a constraint (given in printf form) to the SoS program.
#     This will be the main user-facing way to add constraints.
macro constraint(sos,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        constraint( $(esc(sos)), str )
    end
end

macro partconstraint(sos,key,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        partconstraint( $(esc(sos)), $(esc(key)), str )
    end
end


#  Set a polynomial (given in printf form) as the objective.
macro objective(sos,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        objective( $(esc(sos)), str )
    end
end

#  Add a polynomial (given in printf form) to the objective.
macro partobjective(sos,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        partobjective( $(esc(sos)), str )
    end
end

