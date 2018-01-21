
# Path to rust library relative this file.
libpath = *(@__DIR__,"/../target/debug/")
push!(Libdl.DL_LOAD_PATH,libpath)

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
function safe_string(instring::String, vars::Array{String,1}=VARS)
    # regex from http://www.regular-expressions.info/floatingpoint.html
    for v in vars
        regex = "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
        regex = regex * "$v"
        regex = Regex(regex)
        while typeof(match(regex, instring)) != Void
            m = match(regex, instring).match[1:end-1]
            instring = replace(instring, regex, m * " * $v", 1)
        end
    end
    return instring
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
function init_elist(head::Array{String,1}=HEAD, head_l::Int=HEAD_L, tail::Array{String,1}=TAIL, tail_l::Int=TAIL_L)
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
function parse_elist(elist::String=init_elist(), dict::Dict=DICT)
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
function parse_tree(elist::String=init_elist(), operators=OPERATORS, functions=FUNCTIONS, terminators=TERMINATORS)
    #TODO: Add types to the above!
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
function init_gene(head::Array{String,1}=HEAD, head_l::Int=HEAD_L, tail::Array{String,1}=TAIL, tail_l::Int=TAIL_L, dict::Dict=DICT)
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

"""
    init_chromo(glen, head, head_l, tail, tail_l, dict)

Initialise a Chromosome.

# Examples
```julia-repl
julia> init_chromo(5, ["*"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((2)+(x))+(1))*(((1)+(2))+(2))*(((1)+(2))+(2))*(((1)+(x))+(2))*(((x)+(x))+(2))

```
"""
function init_chromo(glen::Int=GLEN, header_operators::Array{String,1}=HEADER_OPERATORS, head::Array{String,1}=HEAD, head_l::Int=HEAD_L, tail::Array{String,1}=TAIL, tail_l::Int=TAIL_L, dict::Dict=DICT)
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

########################
# Individual functions
########################

#XXX: Rewrite this example
"""
    gen_indi_attributes(clist, de, bc, ival[, flist])

Generates all attributes for an `Individual` given a `clist::Array{Chromosome,1}`.

# Examples
```julia-repl
julia>

```
"""
function gen_indi_attributes(indi_clist::Array{Chromosome,1}, de::Array{String,1}, bc, ival::Array{Tuple{Float64,Float64},1}, flist::Array{String,1}=FLIST, vars=VARS)

    # Prepare for mevac calls
    len = length(vars)
    mevac_vars = map(x -> Base.unsafe_convert(Cstring, Base.cconvert(Cstring, x)), vars)

    # Error (Differential equation)
    ## Set up domain
    dom = map(x -> linspace(x...,10), ival)
    dom = product(dom...)
    dom = unique(dom)

    ## Plug in the expressions in all differential equations and form a sum-of-squares
    # The `de` is given as a list of string-reps.~of the differential equations.
    # Square and add.
    indi_def = parse_expr(de, indi_clist, flist)
    indi_def = *("(" * indi_def[1] * ")^2", map(x-> " + ($x)^2", indi_def[2:end])...)

    ## Clean it up for rust/meval
    # silly silly rust using /wrong/ notation
    indi_def = replace(indi_def, "log", "ln")
    # make raname `Inf` and `NaN` for rust
    indi_def = replace(indi_def, "Inf", "1/0.0")
    indi_def = replace(indi_def, "NaN", "0.0/0.0")
    # remove juxtaposed expressions
    indi_def = safe_string(indi_def)

    # mevac calls
    indi_error = 0
    for pt in dom
        indi_error += ccall(
            (:evalpt,"libmevac"),
            Float64,
            (Cstring,Array{Cstring,1},Array{Float64,1},UInt32),
            indi_def, mevac_vars, [pt...], len)
    end

    # Penalty (Boundary condition)
    # The boundary conditions are given as an array of [string-rep., [point]] so we
    # loop over all the boundary conditions at these points and square and add.
    indi_penalty = 0
    lambda = 100
    for b in bc
        # Set up a domain
        # For PDEs, boundary conditions are not points,
        # so need to be calculated over domain instead.
        dom = Tuple[]
        for i in 1:length(b[2])
            if typeof(b[2][i]) == String
                dom = vcat(dom, tuple(ival[i]...))
            else
                dom = vcat(dom, tuple(b[2][i], b[2][i]))
            end
        end
        dom = map(x -> linspace(x..., 10), dom)
        dom = product(dom...)
        dom = unique(dom)

        # Plug in, and square expression
        indi_bcf = parse_expr(b[1], indi_clist, flist)
        indi_bcf = "(" * indi_bcf * ")^2"
        # Clean it up for rust/meval
        indi_bcf = replace(indi_bcf, "log", "ln")
        indi_bcf = replace(indi_bcf, "Inf", "1/0.0")
        indi_bcf = replace(indi_bcf, "NaN", "0.0/0.0")
        indi_bcf = safe_string(indi_bcf)

        # call mevac
        for pt in dom
            indi_penalty += ccall(
                (:evalpt,"libmevac"),
                Float64,
                (Cstring,Array{Cstring,1},Array{Float64,1},UInt32),
                indi_bcf, mevac_vars, [pt...], len)
        end
    end
    indi_penalty = lambda*indi_penalty

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
function init_indi(de::Array{String,1}, bc, ival::Array{Tuple{Float64,Float64},1}, flist::Array{String,1}=FLIST, glen::Int=GLEN, header_operators::Array{String,1}=HEADER_OPERATORS, head::Array{String,1}=HEAD, head_l::Int=HEAD_L, tail::Array{String,1}=TAIL, tail_l::Int=TAIL_L, dict::Dict=DICT)
    # List of Chromosomes
    indi_clist = map(x -> init_chromo(glen, header_operators, head, head_l, tail, tail_l, dict), 1:length(flist))
    return Individual(indi_clist, gen_indi_attributes(indi_clist, de, bc, ival, flist)...)
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
function reparse_indi(inindi::Individual, flist::Array{String,1}=FLIST)
    indi = deepcopy(inindi)
    indi_clist = indi.clist

    de = indi.de
    bc = indi.bc
    ival = indi.ival

    return Individual(indi_clist, gen_indi_attributes(indi_clist, de, bc, ival, flist)...)
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

#TODO: Rewrite all this documentation.
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

WARN: This function does not handle systems correctly at this point.

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
    m = matchall(regexp, expr)
    for i in m
        if i[1] == 'e'
            sub = thestring
            for var in i[2:end]
                try
                    sub = differentiate(sub, parse(string(var)))
                catch
                    sub = "Inf"
                end
            end
            expr = replace(expr, "<" * i * ">", sub)
        end
    end
    return expr
end

function parse_expr(expr::String, chromos::Array{Chromosome,1}, flist::Array{String,1}=FLIST)
    zipped = zip(flist, chromos)
    for z in zipped
        te = replace(expr, z[1], "e")
        expr = parse_expr(te, z[2])
    end
    return expr
end

function parse_expr(expr::String, indi::Individual, flist::Array{String,1}=FLIST)
    return parse_expr(expr, indi.clist, flist)
end

function parse_expr(exprs::Array{String,1}, chromos::Array{Chromosome,1}, flist::Array{String,1}=FLIST)
    return map(x -> parse_expr(x, chromos, flist), exprs)
end

function parse_expr(exprs::Array{String,1}, indi::Individual, flist::Array{String,1}=FLIST)
    return map(x -> parse_expr(x, indi.clist, flist), exprs)
end
