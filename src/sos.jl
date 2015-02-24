#############################
###### Sums of Squares ######
#############################

#  A system of polynomial equations.
#    We convert inequalities into equations with slack variables (squared).
type SoS
    slackvars::Int  # counter
    vars::Set{Symbol}  # all variables, including slack variables
    constraints::Array{SoSPoly}  # polynomials that are constrained to 0
    objective::SoSPoly
end
SoS() = SoS(0, Set{Symbol}(), [], SoSPoly())

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




#  Parse a polynomial (in Expr form) into a SoSPoly, by traversing
#    a Julia syntax tree.
function parsepoly(sos::SoS, ex::Symbol)
    push!(sos.vars, ex)
    return [[ex => 1] => 1.0] :: SoSPoly
end
function parsepoly{T<:Number}(sos::SoS, ex::T)
    return [(Symbol=>Int64)[] => convert(Float64, ex)] :: SoSPoly
end
function parsepoly(sos::SoS, ex::Expr)
    if(ex.head == :call)
        if(ex.args[1] == :(+))
            return parsepoly(sos,ex.args[2]) + parsepoly(sos,ex.args[3])
        elseif(ex.args[1] == :(*))
            return parsepoly(sos,ex.args[2]) * parsepoly(sos,ex.args[3])
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
    end
    throw(ArgumentError(@sprintf("Expression not understood as polynomial: %s", ex)))
end


#  Add any constraint (given in Expr form) to the SoS program, by
#    massaging it into a polynomial equaling zero.
function addconstraint(sos :: SoS, ex :: Expr)
    if(ex.head == :comparison)
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
    elseif(ex.head == :(=)) # isn't Julia comparison, but we should support
        polyex = :( $(ex.args[1]) - $(ex.args[2]) )
    else
        throw(ArgumentError("Constraint must be an equation or inequality"))
    end
    #@printf("Original expression: %s\n", ex)
    #@printf("Zeroing polynomial expression: %s\n", polyex)
    poly = parsepoly(sos, polyex)
    #@printf("Parsed poly: %s\n", poly)
    push!(sos.constraints, poly)
end

#  Add a polynomial (given in Expr form) to the objective.
function addobjective(sos :: SoS, ex)
    poly = parsepoly(sos, ex)
    sos.objective += poly
end
#  Set a polynomial (given in Expr form) as the objective, overwriting what
#    was there before.
function setobjective(sos :: SoS, ex)
    poly = parsepoly(sos, ex)
    sos.objective = poly
end

#  The following three macros are the main user-facing ways to prepare
#    a SoS program.
#    This section is `esc` hell and was a pain to figure out.

#  Add a constraint (given in printf form) to the SoS program.
#     This will be the main user-facing way to add constraints.
macro constraint(sos,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        addconstraint( $(esc(sos)), parse(str) )
    end
end

#  Add a polynomial (given in printf form) to the objective.
macro partobjective(sos,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        addobjective( $(esc(sos)), parse(str) )
    end
end

#  Set a polynomial (given in printf form) as the objective.
macro partobjective(sos,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        setobjective( $(esc(sos)), parse(str) )
    end
end

