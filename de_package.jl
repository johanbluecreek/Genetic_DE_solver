using Calculus
#XXX: The above cause warnings if running this file with @everywhere, but seems to work...

using Iterators

import Base.show

#XXX: Importing == from Base and rewriting it for comparing Individuals causes a
# fatal error:
#=
include("./de_package.jl")
WARNING: Method definition ==(Any, Any) in module Base at operators.jl:15 overwritten in module Main at /home/user/github/Genetic_DE_solver/de_package.jl:636.
fatal: error thrown and no exception handler available.
Base.MethodError(f=Main.#eq_indi(), args=(:failed, :failed))
rec_backtrace at /home/user/git_clones/julia/src/stackwalk.c:84
record_backtrace at /home/user/git_clones/julia/src/task.c:233
jl_throw at /home/user/git_clones/julia/src/task.c:551
jl_method_error_bare at /home/user/git_clones/julia/src/gf.c:1091
jl_method_error at /home/user/git_clones/julia/src/gf.c:1108
jl_apply_generic at /home/user/git_clones/julia/src/gf.c:1941
task_done_hook at ./task.jl:144
unknown function (ip: 0x7f88b1ad20c2)
jl_call_method_internal at /home/user/git_clones/julia/src/julia_internal.h:210 [inlined]
jl_apply_generic at /home/user/git_clones/julia/src/gf.c:1950
jl_apply at /home/user/git_clones/julia/src/julia.h:1392 [inlined]
finish_task at /home/user/git_clones/julia/src/task.c:215 [inlined]
start_task at /home/user/git_clones/julia/src/task.c:262
unknown function (ip: 0xffffffffffffffff)
=#
# import Base.==

################################################
#
# Hard-coded global variables
#
################################################

# Basis set
functions = ["s", "c", "e", "l", "u"]
operators = ["+", "-", "*", "/"]
digits = vcat(["$i" for i=range(0,10)], ["p"])
# Relevant variables for PDEs (can also be used as unspecified constants)
vars = ["x", "y"]
# Relevant function labels for systems
flist = ["f1", "f2"]

terminators = vcat(digits, vars, vars, vars, vars, vars)

# "z" is an operator that deactivates the following chromosome
header_operators = vcat(operators, ["z"])

head = vcat(functions, operators, vars, digits, vars, vars, vars, vars)
tail = vcat(digits, vars)

dict = Dict(
  "+" => "(<expr>)+(<expr>)",
  "-" => "(<expr>)-(<expr>)",
  "*" => "(<expr>)*(<expr>)",
  "/" => "(<expr>)/(<expr>)",
  "s" => "sin(<expr>)",
  "c" => "cos(<expr>)",
  "l" => "log(<expr>)",
  "e" => "exp(<expr>)",
  "u" => "(<expr>)",
  "p" => "pi"
)

################################################
#
# Types
#
################################################

"""
    Gene

Most basic type holding one mathematica expression.

Represents small/atomic mathematical expressions that should combine into expressions
representing solutions in a `Chromosome`.
"""
type Gene
  #XXX: Should act as old Chromosome
  elist::String         # String representation of the gene
  thestring::String     # String representation of the expression
  tree::String          # Tree-like representation of the clist

  #XXX: See where it is most appropriate to save these.
  head::Array{String,1}             # Save of head list used
  head_l::Int                       # Save of head length used
  tail::Array{String,1}             # Save of tail list used
  tail_l::Int                       # Save of tail length used

  dict::Dict                        # Save of dictionary used
end
show(io::IO, x::Gene) = print(io, x.thestring)

"""
    Chromosome

A type holding several `Gene` types to form a more complex mathematical expression.

Represents a complete solution to a differential equation, possibly in a system. These are
mathematical expressions independent of the differential equation in these sense of that it
does not have a fitness. A collection of `Chromosome` types combine to an `Individual` that
represents the solution to the whole system (note that the system can be a system of one).
"""
type Chromosome
  #XXX: Should act as old Individual (except fitness)
  glist::Array{Gene,1}    # List of Genes

  header::String          # To combine Genes

  thestring::String       # Full string representation of mathematical expression
end
show(io::IO, x::Chromosome) = print(io, x.thestring)

"""
    Individual

Types representing complete solutions to systems of differential equations.
"""
type Individual
  clist::Array{Chromosome,1}  # List of Chromosomes
  def::Function               # Function to evaluate the de

  # Fitness and its components
  error::Float64
  penalty::Float64
  shape::Float64
  fitness::Float64

  # For save
  de::Array{String,1}
  bc::Array{Any,1}
  ival::Array{Tuple{Float64,Float64},1}
end
show(io::IO, x::Individual) = print(io, string(x.fitness) * ": " * join(x.clist, ", "))

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

########################
# "elist"-functions
########################

"""
    init_elist(head, head_l, tail, tail_l)

Initialise an elist for a `Gene`.

Note that `tail_l` should be twice as long as `head_l` (or longer) to ensure that the
mathematical expression truncates. Choose it smaller to live on the wild side.

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
function init_chromo(glen::Int=2, header_operators::Array{String,1}=header_operators, head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
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
 ```julia-repl
julia> chromo = init_chromo(2, ["*"], ["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)", "*" => "(<expr>)*(<expr>)"))
(((1)+(x))+(2))*(((x)+(2))+(1))

julia> expr = "<e>"; parse_expr(expr, chromo)
"(((1)+(x))+(2))*(((x)+(2))+(1))"

julia> expr = "<ex>"; parse_expr(expr, chromo)
"(1 + (2 + x)) + (2 + (1 + x))"

```

## System of "PDEs" (advanced example)
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
  # Error (Differential equation)
  ## Differential equation as a function (sum of squares)
  indi_def = parse_expr(de, indi_clist, flist)
  indi_def = *(map(x-> " + ($x)^2", indi_def)...)
  indi_def = "(" * join(vars, ",") * ") -> " * indi_def
  indi_def = eval(parse(indi_def))
  ## Set up domain
  dom = map(x -> linspace(x...,10), ival)
  dom = product(dom...)
  ## Calculate error
  indi_error = +(map(x -> indi_def(x...), dom)...)

  # Penalty (Boundary condition)
  # XXX
  indi_penalty = 0

  # Shape-error (Higher derivative)
  # XXX
  indi_shape = 0

  # Fitness
  indi_fitness = indi_error + indi_penalty + indi_shape

  # Return
  return indi_def, indi_error, indi_penalty, indi_shape, indi_fitness, de, bc, ival
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
function init_indi(de::Array{String,1}, bc, ival::Array{Tuple{Float64,Float64},1}, flist::Array{String,1}=flist, glen::Int=2, header_operators::Array{String,1}=header_operators, head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
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

########################
# Genetic functions and operators
########################

############
# Population functions
############

"""
    gen_pop(pop_size, de, bc, ival[, flist, glen, header_operators, head, head_l, tail, tail_l, dict])

Generates a population (`Array{Individual,1}`).

# Examples
```julia-repl
julia> de = ["<f1x> + <f2x>", "<f1y> - (<f2y>)"];

julia> bc = [1];

julia> ival = [(0.1,1.0), (0.2,1.1)];

julia> head_l = 2; tail_l = 4;

julia> gen_pop(5, de, bc, ival)
5-element Array{Individual,1}:
 900.0: (3)*(y), (log(3))-(pi)
 3003.6314008246723: (x)/(9), (7)*(sin(y))
 319.0: (x), (x)*(y)
 2500.0: (5)*(y), (6)
 100.0: ((2)/(pi))+(log(cos(8))), (x)

```
"""
function gen_pop(pop_size::Int, de::Array{String,1}, bc, ival::Array{Tuple{Float64,Float64},1}, flist::Array{String,1}=flist, glen::Int=2, header_operators::Array{String,1}=header_operators, head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
  pop = Individual[]
  while length(pop) < pop_size
    indi = init_indi(de, bc, ival, flist, glen, header_operators, head, head_l, tail, tail_l, dict)
    same = false
    for mem in pop
      if eq_indi(indi, mem)
        same = true
      end
    end
    if !same
      push!(pop, indi)
    end
  end
  return pop
end

"""
    sort_pop(pop)

Sorts a `pop::Array{Individual,1}` according to fitness.

This is simply a short-hand for `sort(pop, by=x -> x.fitness)`.

"""
function sort_pop(pop::Array{Individual,1})
  return sort(pop, by=x -> x.fitness)
end

############
# Mutation functions
############

"""
    mut_change(gene::Gene, mrate::Float64)

Changes entry to one of the same type.

`mrate` is the likelyhood that an entry will change, e.g. `rand() <= mrate`
determines change.

# Examples
```julia-repl
julia> operators = ["+"];

julia> terminators = ["1", "2", "x"];

julia> gene = init_gene(["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)"))
((2)+(x))+(1)

julia> mut_change(gene, 1.0)
((x)+(2))+(x)

julia> mut_change(gene, 1.0)
((1)+(1))+(1)

```
"""
function mut_change(gene::Gene, mrate::Float64)
  new_elist = ""
  for part in gene.elist
    if rand() <= mrate
      #XXX: This does not work great with global variables like this...
      if string(part) in functions
        new_elist *= rand(functions)
      elseif string(part) in operators
        new_elist *= rand(operators)
      elseif string(part) in terminators
        new_elist *= rand(terminators)
      else
        println("WARNING: '$part' not found in the lists 'functions', 'operators', nor 'terminators'. Nothing done.")
        new_elist *= string(part)
      end
    else
      new_elist *= string(part)
    end
  end
  return reparse_gene(Gene(new_elist, gene.thestring, gene.tree, gene.head, gene.head_l, gene.tail, gene.tail_l, gene.dict))
end

"""
    mut_random(gene::Gene, mrate::Float64)

Changes entry to a random entry (respecting head/tail).

`mrate` is the likelyhood that an entry will change, e.g. `rand() <= mrate`
determines change.

# Examples
```julia-repl
julia>

```
"""
function mut_random(gene::Gene, mrate::Float64)
  new_elist = ""
  iter = 1
  for part in gene.elist
    if rand() <= mrate
      if iter <= gene.head_l
        new_elist *= rand(gene.head)
      else
        new_elist *= rand(gene.tail)
      end
    else
      new_elist *= string(part)
    end
    iter += 1
  end
  return reparse_gene(Gene(new_elist, gene.thestring, gene.tree, gene.head, gene.head_l, gene.tail, gene.tail_l, gene.dict))
end

"""
    mut_trunc(gene::Gene)

Truncates a branch in a `Gene` expression tree, and pads the end of the `elist`.

# Examples
```julia-repl
julia>

```
"""
function mut_trunc(genein::Gene)
  gene = deepcopy(genein)

  # Truncation only acts on genes with more than one branch
  if length(gene.tree) == 3
    return gene
  end

  branches = find( x -> x == '(', gene.tree)

  # Select a random branch
  b_start = rand(branches)
  # construct the branch
  i = b_start+1
  level = 1
  branch = "("
  while level != 0
    entry = gene.tree[i]
    branch *= string(entry)
    i += 1
    if entry == '('
      level += 1
    elseif entry == ')'
      level -= 1
    end
  end

  # truncating something that is already a terminator is no use
  if length(branch) == 3
    return gene
  end

  # flatten it
  branch = replace(replace(branch,"(",""),")","")

  #XXX: Note that this might cause unwanted replacements
  new_elist = replace(gene.elist, branch, rand(gene.tail))

  # Pad it to cannonical length
  while length(new_elist) < length(gene.elist)
    new_elist *= string(rand(gene.tail))
  end

  return reparse_gene(Gene(new_elist, gene.thestring, gene.tree, gene.head, gene.head_l, gene.tail, gene.tail_l, gene.dict))
end

function mut_trunc(genein::Gene, mrate::Float64)
    return mut_trunc(genein)
end

"""
    mut_grow(gene::Gene)

Grows a branch in a `Gene` expression tree from a terminator, and truncates from the end.

# Examples
```julia-repl
julia>

```
"""
function mut_grow(genein::Gene)
  gene = deepcopy(genein)

  # Find a terminator that is to be replaced by a tree
  termins = find( x -> string(x) in gene.tail, gene.tree)
  pos = rand(termins)

  # Make a new tree
  #XXX: Can be a long loop... (but seems to work ok)
  new_tree = ""
  while length(new_tree) <= 3
    new_tree = init_elist(gene.head, gene.head_l, gene.tail, gene.tail_l)
    new_tree = parse_tree(new_tree)
  end

  # Insert new tree
  new_elist = String[]
  for i in 1:length(gene.tree)
    if i == pos
      push!(new_elist, new_tree)
    else
      push!(new_elist, string(gene.tree[i]))
    end
  end
  new_elist = string(new_elist...)
  new_elist = replace(replace(new_elist,"(",""),")","")

  # complete the elist
  old_elist = replace(replace(gene.tree,"(",""),")","")
  new_elist = replace(gene.elist, old_elist, new_elist)

  # Make cannoncical length
  new_elist = new_elist[1:length(gene.elist)]

  # Make it safe
  new_elist = ""
  for i in 1:length(new_elist)
    char = string(new_elist[i])
    if char in gene.head && !(char in gene.tail) && i > gene.head_l
      nnew_elist *= rand(gene.tail)
    else
      nnew_elist *= string(new_elist[i])
    end
  end
  new_elist = nnew_elist

  # return
  return reparse_gene(Gene(new_elist, gene.thestring, gene.tree, gene.head, gene.head_l, gene.tail, gene.tail_l, gene.dict))
end

function mut_grow(genein::Gene, mrate::Float64)
    return mut_grow(genein)
end

"""
    mut_swap(gene::Gene)

Swaps the arguments for an operator, e.g. "/12" becomes "/21".

# Examples
```julia-repl
julia> functions = String[];

julia> operators = ["+"];

julia> terminators = ["1", "2", "x"];

julia> gene = init_gene(["+"], 2, ["1", "2", "x"], 4, Dict("+" => "(<expr>)+(<expr>)"))
((2)+(2))+(x)

julia> mut_swap(gene)
(x)+((2)+(2))

```
"""
function mut_swap(genein::Gene)
  gene = deepcopy(genein)

  # Find a operator that is to have its arguments swaped
  ops = find( x -> string(x) in operators, gene.tree)
  # return original if there are no operators
  if length(ops) == 0
    return gene
  end
  pos = rand(ops)
  # construct first sub-tree
  fts = pos+1
  open = 1
  c = fts+1
  while open > 0
    if gene.tree[c] == '('
      open += 1
    elseif gene.tree[c] == ')'
      open -= 1
    end
    c += 1
  end
  fte = c-1
  ft = gene.tree[fts:fte]
  # construct second sub-tree
  sts = fte+1
  open = 1
  c = sts+1
  while open > 0
    if gene.tree[c] == '('
      open += 1
    elseif gene.tree[c] == ')'
      open -= 1
    end
    c += 1
  end
  ste = c-1
  st = gene.tree[sts:ste]

  newtree = gene.tree[1:pos]
  newtree *= st
  newtree *= ft
  newtree *= gene.tree[c:end]

  newtree = replace(replace(newtree,"(",""),")","")
  origtree = replace(replace(gene.tree,"(",""),")","")
  #XXX: This is may cause unwanted replacements
  new_elist = replace(gene.elist, origtree, newtree)

  return reparse_gene(Gene(new_elist, gene.thestring, gene.tree, gene.head, gene.head_l, gene.tail, gene.tail_l, gene.dict))
end

function mut_swap(genein::Gene, mrate::Float64)
    return mut_swap(genein)
end

"""
    mutate(input, mrate, mselchance, method)

Mutates `input` using `method`.

# Inputs

The `input` can be either `Chromosome` or `Individual`, where for an individual all
chromosome parts will be treated equally.

# Methods

`["change", "swap", "grow", "trunc", "random"]`

See documentation for `mut_[method]` for a description for each of the methods.

# Examples
```julia-repl
julia> de = ["<f1x> + <f2x>", "<f1y> - (<f2y>)"];

julia> bc = [1];

julia> ival = [(0.1,1.0), (0.2,1.1)];

julia> head_l = 2; tail_l = 4;

julia> indi = init_indi(de, bc, ival)
3700.0: (6)*(y), ((y)+(x))-(y)

julia> mutate(indi)
2500.0: (6)*(y), ((y)+(x))-(x)

julia> mutate(indi)
500.0: (2)*(y), ((y)+(x))-(y)

```
"""
function mutate(inchromo::Chromosome, mrate::Float64=0.6, mselchance::Float64=0.2, method::String="change")
    chromo = deepcopy(inchromo)
    methods = ["change", "swap", "grow", "trunc", "random"]
    if method in methods
        if length(chromo.glist) == 1 #XXX: Not been tested yet!
            ee = "mut_$method"
            ee = eval(parse(ee))
            Gene[ee(chromo.glist[1], mrate)]
            chromo.glist = eval(parse(ee))
            chromo = reparse_chromo(chromo)
            return chromo
        else
            new_glist = Gene[]
            for gene in chromo.glist
                if rand() <= mselchance
                    ee = "mut_$method"
                    ee = eval(parse(ee))
                    ee = ee(gene, mrate)
                    push!(new_glist, ee)
                else
                    push!(new_glist, gene)
                end
            end
            chromo.glist = new_glist
            chromo = reparse_chromo(chromo)
            return chromo
        end
    else
        println("WARNING: No support for method: '$method'. Nothing done.")
        return chromo
    end
end

function mutate(inindi::Individual, mrate::Float64=0.6, mselchance::Float64=0.2, method::String="change")
    indi = deepcopy(inindi)
    new_clist = Chromosome[]
    for chromo in indi.clist
        push!(new_clist, mutate(chromo, mrate, mselchance, method))
    end
    indi.clist = new_clist
    indi = reparse_indi(indi)
    return indi
end

############
# Crossover functions
############

###########################################
# XXX OLD STUFF XXX
###########################################

### GP Operators ###

# Mutate -- Chromosome level

## Deps for mutate

## Actual mutate
#TODO: Finish the other methods
#TODO: Some meta-programming in the below should reduce the code-length significantly e.g. eval(parse("mut_" * method))
#=
function mutate(inindi::Individual, mrate::Float64=0.6, mselchance::Float64=0.2, method::String="change")
  indi = deepcopy(inindi)
  if method == "change"
    # Change entry to entry of same type, see mut_change()
    if length(indi.clist) == 1
      indi.clist = Chromosome[mut_change(indi.clist[1], mrate)]
      indi = reparse_indi(indi)
      return indi
    else
      new_clist = Chromosome[]
      for chromo in indi.clist
        if rand() <= mselchance
          push!(new_clist, mut_change(chromo, mrate))
        else
          push!(new_clist, chromo)
        end
      end
      indi.clist = new_clist
      indi = reparse_indi(indi)
      return indi
    end
  elseif method == "swap"
    # Swaps arguements for an operator
    if length(indi.clist) == 1
      indi.clist = Chromosome[mut_swap(indi.clist[1])]
      indi = reparse_indi(indi)
      return indi
    else
      new_clist = Chromosome[]
      for chromo in indi.clist
        if rand() <= mselchance
          push!(new_clist, mut_swap(chromo))
        else
          push!(new_clist, chromo)
        end
      end
      indi.clist = new_clist
      indi = reparse_indi(indi)
      return indi
    end
  elseif method == "grow"
    # Grows a random tree
    if length(indi.clist) == 1
      indi.clist = Chromosome[mut_grow(indi.clist[1])]
      indi = reparse_indi(indi)
      return indi
    else
      new_clist = Chromosome[]
      for chromo in indi.clist
        if rand() <= mselchance
          push!(new_clist, mut_grow(chromo))
        else
          push!(new_clist, chromo)
        end
      end
      indi.clist = new_clist
      indi = reparse_indi(indi)
      return indi
    end
    return indi
  elseif method == "trunc"
    # Truncates a random tree to a random terminator
    if length(indi.clist) == 1
      indi.clist = Chromosome[mut_trunc(indi.clist[1])]
      indi = reparse_indi(indi)
      return indi
    else
      new_clist = Chromosome[]
      for chromo in indi.clist
        if rand() <= mselchance
          push!(new_clist, mut_trunc(chromo))
        else
          push!(new_clist, chromo)
        end
      end
      indi.clist = new_clist
      indi = reparse_indi(indi)
      return indi
    end
    return indi
  elseif method == "random"
    # Changes entries randomly
    if length(indi.clist) == 1
      indi.clist = Chromosome[mut_random(indi.clist[1], mrate)]
      indi = reparse_indi(indi)
      return indi
    else
      new_clist = Chromosome[]
      for chromo in indi.clist
        if rand() <= mselchance
          push!(new_clist, mut_random(chromo, mrate))
        else
          push!(new_clist, chromo)
        end
      end
      indi.clist = new_clist
      indi = reparse_indi(indi)
      return indi
    end
  else
    println("WARNING: No support for method: '$method'. Nothing done.")
    return indi
  end
end
=#
# Mutate -- Individual level

#TODO: Write a mutate that changes the places of chromosomes in the individual.
#XXX: Random change for now. Implement support for mrate
function mutate_jump(inindi::Individual, mrate::Float64=0.3)
  indi = deepcopy(inindi)
  new_clist = map(x->indi.clist[x], randperm(length(indi.clist)))

  #XXX: The below code-block is repeated in three places, at least! Make it a function!
  thestring = ""
  for i in (length(indi.clist)-1):-1:1
    if string(indi.header[i]) != "z"
      thestring = string(indi.header[i]) * "(" * new_clist[i+1].thestring * ")" * thestring
    end
  end
  thestring = "(" * new_clist[1].thestring * ")" * thestring

  return init_full_indi(
    Individual(new_clist, indi.header, thestring, Inf, Inf, Inf, Inf, "", ["",0], [0.0,0.0], indi.head, indi.head_l, indi.tail, indi.tail_l, indi.dict),
    indi.de, indi.bc, indi.ival
  )
end

# Mutates the header that combines the Chromosomes beloning to an Individual
function mutate_head(inindi::Individual, mrate::Float64=0.3)
  indi = deepcopy(inindi)
  new_header = ""
  for part in indi.header
    if rand() <= mrate
      new_header *= rand(operators)
    else
      new_header *= string(part)
    end
  end
  # generate new thestring
  # TODO: Merge the below with the above loop, to avoid having two for-loops (?)
  new_thestring = ""
  for i in (length(indi.clist)-1):-1:1
    if string(new_header[i]) != "z"
      new_thestring = string(new_header[i]) * "(" * indi.clist[i+1].thestring * ")" * new_thestring
    end
  end
  new_thestring = "(" * indi.clist[1].thestring * ")" * new_thestring

  return init_full_indi(
    Individual(indi.clist, new_header, new_thestring, Inf, Inf, Inf, Inf, "", ["",0], [0.0,0.0], indi.head, indi.head_l, indi.tail, indi.tail_l, indi.dict),
    indi.de, indi.bc, indi.ival
  )
end

# Crossover -- Chromosome level

## Actual crossover
#TODO: Implement the other methods
function crossover(p1in::Individual, p2in::Individual, cselchance::Float64=1.0, method::String="one-point")
  p1, p2 = deepcopy(p1in), deepcopy(p2in)
  if method == "one-point"
    # Breaks a Chromosome clist at a random point. Causes chaotic behaviour.
    for i in 1:length(p1.clist)
      if rand() <= cselchance
        new_cclist1 = ""
        new_cclist2 = ""
        # Use crossover in active genes in head so that it actually gives change
        l1 = length(replace(replace(p1.clist[i].tree,"(",""),")",""))
        l2 = length(replace(replace(p2.clist[i].tree,"(",""),")",""))
        rnd = 1
        if l1 <= l2
          rnd = rand(1:l1)
        elseif l2 < l1
          rnd = rand(1:l2)
        else
          # If the above fails, change deactive head
          #XXX: When does this happen?
          rnd = rand(2:head_l)
        end
        new_cclist1 = p1.clist[i].clist[1:rnd] * p2.clist[i].clist[rnd+1:end]
        new_cclist2 = p2.clist[i].clist[1:rnd] * p1.clist[i].clist[rnd+1:end]
        p1.clist[i].clist = new_cclist1
        p2.clist[i].clist = new_cclist2
      end
    end
    return reparse_indi(p1), reparse_indi(p2)
  elseif method == "random-two-point"
    # Breaks a Chromosome clist at two random points. Causes chaotic behaviour.
    return p1, p2
  elseif method == "safe-two-point"
    # Switches sub-branches between individuals
#=
    for i in 1:length(p1.clist)
      if rand() <= cselchance
        br1 = find( x -> x == '(', chromo.tree)
    # Substitutes a branch.
    branches = find( x -> x == '(', chromo.tree)

    # Select a random branch
    b_start = rand(branches)
    return p1, p2
=#
  elseif method == "random-crossing"
    # Allow genes to jump between positions in the Individuals clist
    #TODO: Loop this to make sure that the returned parents are new
    cclist1 = Chromosome[]
    cclist2 = Chromosome[]
    for i in 1:length(p1.clist)
      if rand() <= cselchance
        push!(cclist1, p2.clist[i])
        push!(cclist2, p1.clist[i])
      else
        push!(cclist1, p1.clist[i])
        push!(cclist2, p2.clist[i])
      end
    end
    p1.clist = cclist1
    p2.clist = cclist2
    return reparse_indi(p1), reparse_indi(p2)
  else
    println("WARNING: No support for method: '$method'. Nothing done.")
    return p1, p2
  end
end

# Crossover -- Individual level

#TODO: Implement this
function crossover_head(p1in::Individual, p2in::Individual, method::String="one-point")
  p1, p2 = deepcopy(p1in), deepcopy(p2in)
  if method == "one-point"
    return p1, p2
  elseif method == "two-point"
    return p1, p2
  else
    println("WARNING: No support for method: '$method'. Nothing done.")
    return p1, p2
  end
end

# Parent-selection

# Function to return two parents for crossover
#TODO: Implement other methods
function p_select(pop::Array{Individual, 1}, method::String="tournament")
  if method == "random"
    p1 = rand(pop)
    p2 = rand(pop)
    while eq_indi(p1, p2)
      p2 = rand(pop)
    end
    return p1, p2
  elseif method == "tournament"
    pop_size = length(pop)
    p1 = rand(pop,rand(2:Int(round(sqrt(pop_size)/2))))
    p1 = sort_pop(p1)[1]
    p2 = rand(pop,rand(2:Int(round(sqrt(pop_size)/2))))
    p2 = sort_pop(p2)[1]
    while eq_indi(p1,p2)
      p1 = rand(pop,rand(2:Int(round(sqrt(pop_size)/2))))
      p1 = sort_pop(p1)[1]
      p2 = rand(pop,rand(2:Int(round(sqrt(pop_size)/2))))
      p2 = sort_pop(p2)[1]
    end
    return p1, p2
  end
end

# Generates a population to use for crossover
function gen_p_pop(pop::Array{Individual, 1}, p_pop_halfsize::Int, method::String="tournament")
  p_pop = Tuple{Individual,Individual}[]
  for i in 1:p_pop_halfsize
    push!(p_pop, p_select(pop, method))
  end
  return p_pop
end


# Mutation selection

# Selects Individuals for mutation
#TODO: Implement other methods
function m_select(pop::Array{Individual, 1}, method::String="random")
  if method == "random"
    p1 = rand(pop)
    return p1
  elseif method == "tournament"
    pop_size = length(pop)
    p1 = rand(pop,rand(2:Int(round(sqrt(pop_size)/2))))
    p1 = sort_pop(p1)[1]
    return p1
  end
end

# Generate a population for mutation
#TODO: Implement a uniq option here
function gen_m_pop(pop::Array{Individual, 1}, m_pop_size::Int, method::String="random")
  m_pop = Individual[]
  for i in 1:m_pop_size
    push!(m_pop, m_select(pop, method))
  end
  return m_pop
end
