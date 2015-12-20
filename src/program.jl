#################################
###### Polynomial programs ######
#################################

#  A system of polynomial equations.
#  We convert inequalities into equations with slack variables (squared).
type Program
    slackvars :: Int  # counter
    vars :: Set{Symbol}  # all variables, including slack variables
    constraints :: Array{SoSPoly}  # polynomials that are constrained to 0
    namedconstraints :: Dict{Any,SoSPoly}
    objective :: SoSPoly
    maximize :: Bool
    symmetries :: Array{Dict{Symbol,Symbol}}
end
Program(; maximize=true, minimize=false) = begin
    mx = (maximize && !minimize) # default to maximize; switch to minimize if there's any sign to do so
    Program(0, Set{Symbol}(), [], Dict{Any,SoSPoly}(), SoSPoly(), mx, [])
end

#  Generate a new slack variable (a Symbol)
slackvar(prog :: Program) = symbol( @sprintf("slack%d", prog.slackvars += 1) )

# convenience
monoms(prog :: Program, d :: Int64) = monoms([prog.vars...], d)

typealias Maybe{T} Union{T,Void}


#  Parse a polynomial (in Expr form) into a SoSPoly, by traversing
#    a Julia syntax tree.

parsepoly{T<:Number}(prog :: Maybe{Program}, ex :: T) = Dict( SoSMonom() => convert(Float64, ex) ) :: SoSPoly

function parsepoly(prog :: Maybe{Program}, ex::Symbol)
    prog == nothing || push!(prog.vars, ex)
    Dict(Dict(ex => 1) => 1.0) :: SoSPoly
end


function parsepoly(prog :: Maybe{Program}, ex :: Expr)
    if ex.head == :call
        if(ex.args[1] == :(+))
            ret = parsepoly(prog, ex.args[2])
            for i in 3:length(ex.args)
                ret += parsepoly(prog, ex.args[i])
            end
            return ret
        elseif(ex.args[1] == :(*))
            ret = parsepoly(prog, ex.args[2])
            for i in 3:length(ex.args)
                ret *= parsepoly(prog, ex.args[i])
            end
            return ret
        elseif(ex.args[1] == :(-))
            if(length(ex.args) == 2) # unary minus
                return parsepoly(prog, :( (-1) * $(ex.args[2]) ))
            elseif(length(ex.args) == 3) # binary minus
                return parsepoly(prog, :( $(ex.args[2]) + (-1) * $(ex.args[3]) ))
            end
        elseif(ex.args[1] == :(^))
            if(! isa(ex.args[3],Int64) || ex.args[3] < 0)
                throw(ArgumentError("Bad exponent: $(args[3]). Exponents must be non-negative integers"))
            elseif(ex.args[3] == 0)
                return parsepoly(prog, 1)
            end
            operand = parsepoly(prog, ex.args[2])
            return prod([operand for i in 1:ex.args[3]])
        end
    elseif ex.head == :block # this happens on the RHS of single-equals
        if length(ex.args) > 2
            throw(ArgumentError(@sprintf("Multi-line block in expression: %s", ex)))
        end
        return parsepoly(prog, ex.args[2]) # assume it's a 1-line block
    end
    throw(ArgumentError("Expression not understood as polynomial: $ex"))
end


#  Massage any constraint (given in Expr form) into a polynomial equaling zero.
function relationToPolynomial(prog :: Program, ex :: Expr)
    if ex.head == :comparison
        if(ex.args[2] == :(>=) || ex.args[2] == :(≥))
            polyex = :( $(ex.args[1]) - $(ex.args[3]) - $(slackvar(prog))^2 )
        elseif(ex.args[2] == :(==))
            polyex = :( $(ex.args[1]) - $(ex.args[3]) )
        elseif(ex.args[2] == :(<=) || ex.args[2] == :(≤))
            polyex = :( $(ex.args[3]) - $(ex.args[1]) - $(slackvar(prog))^2 )
        else
            throw(ArgumentError(
                "Comparison operator must be >=, ≥, ==, =, <=, or ≤"))
        end
    elseif ex.head == :(=) # isn't Julia comparison, but we should support
        polyex = :( $(ex.args[1]) - $(ex.args[2]) )
    else
        throw(ArgumentError("Constraint must be an equation or inequality"))
    end

    parsepoly(prog, polyex)
end




#  Perturb the objective by Gaussian noise with standard deviation `tol`.
#  Each moment up to degree `deg`, in the variables known to `prog`, receives
#  noise.
function perturb_objective(prog,tol,deg)
    for md in 1:deg
        for mon in monoms(prog,md)
            poly = Dict( mon => tol * randn() ) :: SoSPoly
            prog.objective += poly
        end
    end
end

#  Perturb the objective by Gaussian noise with standard deviation `tol`.
#  A random choice of `s` terms (with replacement), from the moments up to
#  degree `deg`, in the variables known to `prog`, receives noise.
function perturb_objective_sparse(prog,tol,s,deg)
    mon = [ monoms(prog,i) for i in 0:deg ]
    mon = vcat(mon...)
    for i in 1:s
        m = rand(mon)
        prog.objective += Dict( m => tol*randn() )
    end
end


#  Add a constraint to the program.
function constraint(prog,str)
    poly = relationToPolynomial(prog, parse(str))
    push!(prog.constraints, poly)
end

#  Add a polynomial to a constraint indexed by `key`.
function partconstraint(prog,key,str)
    poly = parsepoly(prog, parse(str))

    if haskey(prog.namedconstraints, key)
        prev = prog.namedconstraints[key]
        addpoly!(prev,poly)
    else
        push!(prog.constraints, poly)
        prog.namedconstraints[key] = poly
    end
end

#  Set a polynomial as the objective.
objective(prog,str) = prog.objective = parsepoly(prog, parse(str))

#  Add a polynomial to the objective.
partobjective(prog,str) = prog.objective += parsepoly(prog, parse(str))


#  `printf`-style macros for constraints and objectives.
#  This section is `esc` hell and was a pain to figure out.
#  Since we're calling a macro from a macro, we have three different
#  scopes in play at once, and we need to make sure everything gets evaluated
#  in the right scope.


#  Add a constraint (given in printf form) to the program.
macro constraint(prog,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        constraint( $(esc(prog)), str )
    end
end

#  Add a polynomial (given in printf form) to a constraint indexed by `key`.
macro partconstraint(prog,key,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        partconstraint( $(esc(prog)), $(esc(key)), str )
    end
end


#  Set a polynomial (given in printf form) as the objective.
macro objective(prog,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        objective( $(esc(prog)), str )
    end
end

#  Add a polynomial (given in printf form) to the objective.
macro partobjective(prog,args...)
    escargs = [esc(x) for x in args[2:end]]
    quote
        str = @sprintf($(args[1]),$(escargs...))
        partobjective( $(esc(prog)), str )
    end
end

