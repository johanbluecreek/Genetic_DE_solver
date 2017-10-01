
################################################
#
# Functions
#
################################################

########################
# Generic functions
########################

"""
    safe_string(instring)

Makes multiplication explicit in strings of mathematical expressions.

# Examples
```julia-repl
julia> x = 3;

julia> eval(parse("2x"))
6

julia> y = safe_string("2x")
"2 * x"

julia> eval(parse(y))
6

```
"""
function safe_string(instring::String)
    # regex from http://www.regular-expressions.info/floatingpoint.html
    regex = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?x"
    while typeof(match(regex, instring)) != Void
        m = match(regex, instring).match[1:end-1]
        instring = replace(instring, regex, m * " * x", 1)
    end
    return instring
end

"""
    protect(instring[, m])

Protects a string from variable substitution.
"""
function protect(instring, m=true)
    retstring = instring
    if m
        retstring = replace(retstring, "exp", "EXP")
        retstring = replace(retstring, "log", "LOG")
        retstring = replace(retstring, "sin", "SIN")
        retstring = replace(retstring, "cos", "COS")
    else
        retstring = replace(retstring, "EXP", "exp")
        retstring = replace(retstring, "LOG", "log")
        retstring = replace(retstring, "SIN", "sin")
        retstring = replace(retstring, "COS", "cos")
    end
    return retstring
end

function unprotect(instring)
    protect(instring, false)
end

########################
# "elist"-functions
########################

"""
    init_elist(head, head_l, tail, tail_l)

Initialise an elist for a `Gene`.

Note that `tail_l` should be twice as long as `head_l` (or longer) to ensure that the
mathematical expression truncates. Choose it smaller to live on the wild side (and
experience crashes, eventually).

# Examples
```julia-repl
julia> head = ["s", "c", "e", "l", "u", "+", "-", "*", "/", "x"];

julia> tail = ["1", "2", "x"];

julia> init_elist(head, 5, tail, 10)
"//xlu21221xxxx2"

```
"""
function init_elist(head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l)
    elist = String[]

    for i in 1:head_l
        push!(elist, rand(head))
    end

    for i in 1:tail_l
        push!(elist, rand(tail))
    end

    elist = join(elist)

    return elist
end

"""
    parse_elist(elist, dict)

Parse an elist to a mathematical expression.

# Examples
```julia-repl
julia> elist = "+21";

julia> dict = Dict("+" => "(<expr>)+(<expr>)");

julia> parse_elist(elist, dict)
"(2)+(1)"

```
"""
function parse_elist(elist::String=init_elist(), dict::Dict=dict)
    thestring = "<expr>"
    for i in elist
        thestring = replace(thestring, "<expr>", get(dict, "$i", "$i"), 1)
    end
    return thestring
end

"""
    parse_tree(elist)

Parse an `elist` to a tree-like expression (open parentasis representing new branch).

# Examples
```julia-repl
julia> elist = "+-3*214";

julia> parse_tree(elist)
"(+(-(3)(*(2)(1)))(4))"

```
"""
function parse_tree(elist::String=init_elist())
    tree = "("
    open = 1
    to_open = 0
    level = Int[]

    iterator = 1

    while open != 0
        if string(elist[iterator]) in operators
            tree = tree * string(elist[iterator]) * "("
            to_open += 1
            push!(level, open)
            open += 1
        end
        if string(elist[iterator]) in functions
            tree = tree * string(elist[iterator]) * "("
            open += 1
        end
        if string(elist[iterator]) in terminators
            tree = tree * string(elist[iterator]) * ")"
            open -= 1
            if to_open > 0
                while level[end] != open
                    tree = tree * ")"
                    open -= 1
                end
                tree = tree * "("
                open += 1
                to_open -= 1
                if to_open != 0
                    try
                        level = level[1:end-1]
                    catch
                        println("This caused an error:")
                        println("level: $level")
                        println("tree: $tree")
                        println("open: $open; to_open: $to_open")
                    end
                else
                    level = Int[]
                end
            else
                while open != 0
                    tree = tree * ")"
                    open -= 1
                end
            end
        end
        iterator += 1
        if iterator == 100
            println("Error: parse_tree: endless loop?")
        end
    end
    return tree
end

########################
# Gene functions
########################

"""
    init_gene(head, head_l, tail, tail_l, dict)

Initialise a Gene.

# Examples
```julia-repl
julia> init_gene(["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)"))
((x)+(2))+(1)

```
"""
function init_gene(head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
    elist = init_elist(head, head_l, tail, tail_l)
    thestring = parse_elist(elist, dict)
    tree = parse_tree(elist)
    return Gene(elist, thestring, tree, head, head_l, tail, tail_l, dict)
end

"""
    reparse_gene(gene)

Re-parse a `Gene` who has had its `elist` changed.

Returns a `Gene` in which the `tree` and `thestring` are rebuilt.

# Examples
```julia-repl
julia> gene = init_gene(["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)"))
((1)+(x))+(2)

julia> gene.elist
"++1x21"

julia> gene.elist = "++x1x2"
"++x1x2"

julia> new_gene = reparse_gene(gene)
((x)+(1))+(x)

```
"""
function reparse_gene(genein::Gene)
    gene = deepcopy(genein)
    gene.thestring = parse_elist(gene.elist, gene.dict)
    gene.tree = parse_tree(gene.elist)
    return gene
end

########################
# Chromosome functions
########################

#XXX: Replaces: init_flat_indi
"""
    init_chromo(glen, head, head_l, tail, tail_l, dict)

Initialise a Chromosome.

# Examples
```julia-repl
julia> init_chromo(5, ["*"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((2)+(x))+(1))*(((1)+(2))+(2))*(((1)+(2))+(2))*(((1)+(x))+(2))*(((x)+(x))+(2))

```
"""
function init_chromo(glen::Int=5, header_operators::Array{String,1}=header_operators, head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
    glist = Gene[]
    thestring = ""
    header = ""
    if glen > 1
        header = *(map(x->rand(header_operators), 1:(glen-1))...)
    end
    glist = map(x->init_gene(head, head_l, tail, tail_l, dict), 1:glen)
    # The "z" operator is resoved from the end.
    for i in (glen-1):-1:1
        if string(header[i]) != "z"
            thestring = string(header[i]) * "(" * glist[i+1].thestring * ")" * thestring
        end
    end
    thestring = "(" * glist[1].thestring * ")" * thestring
    return Chromosome(glist, header, thestring)
end

"""
    reparse_chromo(chromo)

Re-parse a `Chromosome` who has had its `glist` or `header` changed.


# Examples
```julia-repl
julia> chromo1 = init_chromo(5, ["*"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((1)+(x))+(1))*(((2)+(1))+(1))*(((x)+(1))+(2))*(((2)+(x))+(x))*(((1)+(2))+(2))

julia> chromo2 = init_chromo(5, ["/"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "/" => "(<expr>)/(<expr>)"))
(((1)+(1))+(2))/(((x)+(x))+(x))/(((x)+(x))+(2))/(((2)+(x))+(2))/(((1)+(x))+(2))

julia> chromo1.glist = chromo2.glist;

julia> reparse_chromo(chromo1)
(((1)+(1))+(2))*(((x)+(x))+(x))*(((x)+(x))+(2))*(((2)+(x))+(2))*(((1)+(x))+(2))

```
"""
function reparse_chromo(inchromo::Chromosome)
    chromo = deepcopy(inchromo)
    glist = chromo.glist
    glen = length(glist)
    header = chromo.header
    thestring = ""
    for i in (glen-1):-1:1
        if string(header[i]) != "z"
            thestring = string(header[i]) * "(" * glist[i+1].thestring * ")" * thestring
        end
    end
    thestring = "(" * glist[1].thestring * ")" * thestring
    return Chromosome(glist, header, thestring)
end

"""
    eq_chromo(chromo1, chromo2)

Checks if two `Chromosome`s are equal.

Development comment: This should not be seen as a strict equality, but a
similarity.

# Examples
```julia-repl
julia> chromo1 = init_chromo(5, ["*"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((2)+(2))+(2))*(((2)+(x))+(x))*(((2)+(x))+(2))*(((2)+(x))+(2))*(((x)+(1))+(x))

julia> chromo2 = init_chromo(5, ["/"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "/" => "(<expr>)/(<expr>)"))
(((2)+(x))+(2))/(((1)+(2))+(1))/(((x)+(1))+(2))/(((2)+(x))+(x))/(((1)+(x))+(x))

julia> eq_chromo(chromo1, chromo2)
false

julia> eq_chromo(chromo1, chromo1)
true

```
"""
function eq_chromo(chromo1::Chromosome, chromo2::Chromosome)
    # Same expression
    if chromo1.thestring == chromo2.thestring
        return true
    end
    # Ignoring header
    eq = true
    for i in 1:length(chromo1.glist)
        if !(chromo1.glist[i].tree == chromo2.glist[i].tree)
            eq = eq && false
        end
    end
    return eq
end

"""
    parse_expr(expr[s], chromo[s][, flist])

Evaluates a string representation of a mathematical expression w.r.t. the `chromo`.

# Examples

## Single "ODE" (simple example)

The default expression for a function is `e`, so for an ODE/expression with only one
function to be evaluated

 ```julia-repl
julia> chromo = init_chromo(2, ["*"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((1)+(x))+(2))*(((x)+(2))+(1))

julia> expr = "<e>"; parse_expr(expr, chromo)
"(((1)+(x))+(2))*(((x)+(2))+(1))"

julia> expr = "<ex>"; parse_expr(expr, chromo)
"(1 + (2 + x)) + (2 + (1 + x))"

```

## System of "PDEs" (advanced example)

When `flist` is provided to `parse_expr()`, the default is overwritten and the given
function names are used instead

```julia-repl
julia> expr1 = "<f1x> + <f2x>";

julia> expr2 = "<f1y> - (<f2y>)";

julia> chromo1 = init_chromo(2, ["*"], ["+"], 2, ["1", "2", "x", "y"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((y)+(y))+(x))*(((2)+(x))+(1))

julia> chromo2 = init_chromo(2, ["*"], ["+"], 2, ["1", "2", "x", "y"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((y)+(x))+(1))*(((1)+(y))+(1))

julia> exprs = [expr1, expr2];

julia> chromos = [chromo1, chromo2];

julia> parse_expr(exprs, chromos, ["f1", "f2"])
2-element Array{String,1}:
 "(1 + (2 + x)) + ((y + y) + x) + 1 + (1 + y)"
 "2 * (1 + (2 + x)) - ((1 + (1 + y)) + (1 + (y + x)))"

```
"""
function parse_expr(expr::String, chromo::Chromosome)
    thestring = chromo.thestring
    regexp = r"(?<=(\<)).*?(?=\>)"
    m = match(regexp, expr)
    #while typeof(m) != Void
    while (try m.match[1] == 'e' catch typeof(m) != Void; end)
    name = m.match
    sub = thestring
    for var in name[2:end]
        #XXX: You will be able to get an infinite loop here unless you use
        # Calculus with https://github.com/johnmyleswhite/Calculus.jl/pull/107
        try
            sub = differentiate(sub, parse(string(var)))
        catch
            sub = "Inf"
        end
    end
    expr = replace(expr, "<" * name * ">", sub)
    m = match(regexp, expr)
    end
    return expr
end

function parse_expr(expr::String, chromos::Array{Chromosome,1}, flist::Array{String,1}=flist)
  return reduce(
    (x,y) -> parse_expr(replace(x, y[1], "e"), y[2]),
    expr,
    zip(flist, chromos)
  )
end

function parse_expr(expr::String, indi::Individual, flist::Array{String,1}=flist)
    return parse_expr(expr, indi.clist, flist)
end

function parse_expr(exprs::Array{String,1}, chromos::Array{Chromosome,1}, flist::Array{String,1}=flist)
    return map(x -> parse_expr(x, chromos, flist), exprs)
end

function parse_expr(exprs::Array{String,1}, indi::Individual, flist::Array{String,1}=flist)
    return parse_expr(x, indi.clist, flist)
end

########################
# Individual functions
########################

#XXX: Add example below
"""
    gen_indi(clist, de, bc, ival[, length, flist, header_operators, head, head_l, tail, tail_l, dict])

Generates all attributes for an `Individual` given a `clist::Array{Chromosome,1}`.

# Examples
```julia-repl
julia>

```
"""
function gen_indi(indi_clist::Array{Chromosome,1}, de::Array{String,1}, bc, ival::Array{Tuple{Float64,Float64},1}, flist::Array{String,1}=flist)
    #XXX: A lot of time is spent here. Make sure it is optimized as far as possible.

    # Error (Differential equation)
    ## Set up domain
    dom = map(x -> linspace(x...,10), ival)
    dom = product(dom...)

    ## Differential equation as a function (sum of squares)

    exprs = parse_expr(de, indi_clist, flist)
    nexprs = Float64[]

    for e in 1:length(exprs)

        expr = exprs[e]
        expr = protect(expr)

        repall = expr
        for i in 1:length(vars)
            repall = Expr(:call, :replace, repall, vars[i], :(mapvar[$i]))
        end
        repall = Expr(:->, :mapvar, repall)
        repall = eval(repall)
        expr = map(a -> Base.invokelatest(repall, a), tuple(dom...))
        expr = map(unprotect, expr)
        try
            expr = map(x -> eval(parse(x)), expr)
            expr = map(x -> x^2, expr)
            expr = sum(expr)
        catch
            expr = Inf
        end

        push!(nexprs, expr)

    end

    indi_error = sum(nexprs)


    #=
    indi_error = map(a -> Base.invokelatest(repall, a), tuple(dom...))
    println(indi_error)
    =#

    #XXX Get the expression on string form
    #XXX protect the functions names: exp -> EXP replacements for instance
    #XXX have a recursive for loop, that eats all in var...
    #=
        function repall(vars)
           if length(vars) > 1
               return Expr(:call, :replace, :expr, parse(vars[1]), :(x[1]))
           else
               return Expr(:call, :replace, :expr, :(vars[1]), :(x[1]))
            ...
    XXX: Figure it our clever boy

    This should be it:

    map( x -> replace(replace x[1]) x[2], dom) build replace with recursion.

    =#

    #=
    expr = parse_expr(de, indi_clist, flist)
    expr = *(map(x-> " + ($x)^2", expr)...)
    expr = parse(expr)

    erfunc = Expr(:->, Expr(:tuple, (map(parse, vars))...), expr)
    erfunc = eval(erfunc)

    indi_error = Inf
    try
        indi_error = +(map(x -> Base.invokelatest(erfunc, x...), dom)...)
    catch
        indi_error = Inf
    end

    erfunc = nothing
    =#

    #=
    indi_def = parse_expr(de, indi_clist, flist)
    indi_def = *(map(x-> " + ($x)^2", indi_def)...)

    # Define `defunc`
    f_body = parse(indi_def)
    f_call = Expr(:call,:defunc,map(parse, vars)...)
    # eval(Expr(:function,f_call,f_body))
    eval(Expr(:function,f_call,f_body))

    ## Calculate error (try/catch to avoid DomainError)
    indi_error = Inf
    try
        indi_error = +(map(x -> Base.invokelatest(defunc, x...), dom)...)
    catch
        indi_error = Inf
    end
    # XXX: Would be nice to be able to avoid `Base.invokelatest` here.
    =#

    # Penalty (Boundary condition)
    indi_penalty = 0
    #=
    lambda = 100
    for b in bc
        # define bcfunc
        indi_bcf = parse_expr(b[1], indi_clist, flist)
        f_body = parse(indi_bcf)
        f_call = Expr(:call,:bcfunc,map(parse, vars)...)
        eval(Expr(:function,f_call,f_body))
        try
            indi_penalty += (Base.invokelatest(bcfunc, map( x -> x[2] , b[2])...))^2
        catch
            indi_penalty += Inf
        end
    end
    indi_penalty = lambda*indi_penalty
    =#

    # Shape-error (Higher derivative)
    # XXX
    indi_shape = 0

    # Fitness
    indi_fitness = indi_error + indi_penalty + indi_shape

    # Return
    return indi_error, indi_penalty, indi_shape, indi_fitness, de, bc, ival
end



"""
    init_indi(de, bc, ival[, flist, header_operators, head, head_l, tail, tail_l, dict])

Initialise a random individual. Short-hand that uses `gen_indi()`.

# Examples
```julia-repl
julia> de = ["<f1x> + <f2x>", "<f1y> - (<f2y>)"];

julia> bc = [1];

julia> ival = [(0.1,1.0), (0.2,1.1)];

julia> head_l = 2; tail_l = 4;

julia> init_indi(de, bc, ival)
200.0: (x), (y)

```
"""
function init_indi(de::Array{String,1}, bc, ival::Array{Tuple{Float64,Float64},1}, flist::Array{String,1}=flist, glen::Int=5, header_operators::Array{String,1}=header_operators, head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
    # List of Chromosomes
    indi_clist = map(x -> init_chromo(glen, header_operators, head, head_l, tail, tail_l, dict), 1:length(flist))
    return Individual(indi_clist, gen_indi(indi_clist, de, bc, ival, flist)...)
end

"""
    reparse_indi(indi[, flist])

Re-parse an `Individual` if its `clist` have changed. Short-hand for `gen_indi()`.

# Examples
```julia-repl
julia> de = ["<f1x> + <f2x>", "<f1y> - (<f2y>)"];

julia> bc = [1];

julia> ival = [(0.1,1.0), (0.2,1.1)];

julia> head_l = 2; tail_l = 4;

julia> indi1 = init_indi(de, bc, ival)
200.0: (x)+(y), (y)*(2)

julia> indi2 = init_indi(de, bc, ival)
0.0: (0)-(log(log(4))), ((2)*(0))-(7)

julia> indi1.clist[1] = indi2.clist[1]
(0)-(log(log(4)))

julia> reparse_indi(indi1)
400.0: (0)-(log(log(4))), (y)*(2)
```
"""
function reparse_indi(inindi::Individual, flist::Array{String,1}=flist)
    indi = deepcopy(inindi)
    indi_clist = indi.clist

    de = indi.de
    bc = indi.bc
    ival = indi.ival

    return Individual(indi_clist, gen_indi(indi_clist, de, bc, ival, flist)...)
end

"""
    eq_indi(indi1, indi2)

Determines if two `Individual`s are equal.

This is a short-hand running `eq_chromo()` on all `Chromosome`s that are a part of
the `Individual`.

# Examples
```julia-repl
julia> de = ["<f1x> + <f2x>", "<f1y> - (<f2y>)"];

julia> bc = [1];

julia> ival = [(0.1,1.0), (0.2,1.1)];

julia> head_l = 2; tail_l = 4;

julia> indi1 = init_indi(de, bc, ival)
7786.7732432014: (y)*((cos(7))-(9)), (pi)*(x)

julia> indi2 = init_indi(de, bc, ival)
998383.0732874984: ((x))/(y), (exp(pi))/((y)*(2))

julia> eq_indi(indi1, indi2)
false

julia> eq_indi(indi1, indi1)
true

```
"""
function eq_indi(indi1::Individual, indi2::Individual)
    eq = true
    for i in 1:length(indi1.clist)
        if !eq_chromo(indi1.clist[i], indi2.clist[i])
            eq = eq && false
        end
    end
    return eq
end

#XXX: See comment on importing == from Base.
# ==(indi1, indi2) = eq_indi(indi1, indi2)
