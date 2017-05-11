using Calculus
#XXX: The above cause warnings if running this file with @everywhere, but seems to work...

functions = ["s", "c", "e", "l", "u"]
operators = ["+", "-", "*", "/"]
digits = vcat(["$i" for i=range(0,10)], ["p"])
vars = ["x"]

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



# Type for an individual Chromosome
type Chromosome
  clist::String                     # String representation of the chromosome
  thestring::String                 # String rep of expression
  tree::String                      # Tree-like representation of the chromosome

  head::Array{String,1}             # Save of head list used
  head_l::Int                       # Save of head length used
  tail::Array{String,1}             # Save of tail list used
  tail_l::Int                       # Save of tail length used

  dict::Dict                        # Save of dictionary used
end

# Function that constructs clist from head and tail
function init_clist(head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l)
  clist = String[]

  for i in 1:head_l
    push!(clist, rand(head))
  end

  for i in 1:tail_l
    push!(clist, rand(tail))
  end

  clist = join(clist)

  return clist
end

# Function that parse the clist into a string expression
function parse_clist(clist::String=init_clist(), dict::Dict=dict)
  thestring = "<expr>"
  for i in clist
    thestring = replace(thestring, "<expr>", get(dict, "$i", "$i"), 1)
  end
  return thestring
end

# Function that parse the clist into a tree expression
function parse_tree(clist::String=init_clist())
  tree = "("
  open = 1
  to_open = 0
  level = Int[]

  iterator = 1

  while open != 0
    if string(clist[iterator]) in operators
      tree = tree * string(clist[iterator]) * "("
      to_open += 1
      push!(level, open)
      open += 1
    end
    if string(clist[iterator]) in functions
      tree = tree * string(clist[iterator]) * "("
      open += 1
    end
    if string(clist[iterator]) in terminators
      tree = tree * string(clist[iterator]) * ")"
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

# Initialise a chromosome
function init_chromo(head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
  clist = init_clist(head, head_l, tail, tail_l)
  thestring = parse_clist(clist, dict)
  tree = parse_tree(clist)
  return Chromosome(clist, thestring, tree, head, head_l, tail, tail_l, dict)
end

# type for an individual
type Individual
  clist::Array{Chromosome,1}  # List of Chromosomes

  # To combine Chromosomes
  header::String

  # Full expression for individual
  thestring::String

  # Fitness and its components
  fitness::Float64
  error::Float64
  penalty::Float64
  shape_error::Float64

  # For save
  de::String
  bc::Array{Any,1}
  ival::Array{Float64,1}

  head::Array{String,1}             # Save of head list used
  head_l::Int                       # Save of head length used
  tail::Array{String,1}             # Save of tail list used
  tail_l::Int                       # Save of tail length used

  dict::Dict                        # Save of dictionary used
end

# Function to initate a flat, e.g. diff-eq independent, Individual
function init_flat_indi(length::Int=1, head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
  clist = Chromosome[]
  thestring = ""
  header = *(map(x->rand(header_operators), 1:(length-1))...)
  clist = map(x->init_chromo(head, head_l, tail, tail_l, dict), 1:length)
  # The "z" operator is resoved from the end.
  for i in (length-1):-1:1
    if string(header[i]) != "z"
      thestring = string(header[i]) * "(" * clist[i+1].thestring * ")" * thestring
    end
  end
  thestring = "(" * clist[1].thestring * ")" * thestring
  return Individual(clist, header, thestring, Inf, Inf, Inf, Inf, "", ["",0], [0.0,0.0], head, head_l, tail, tail_l, dict)
end

# Function for printing an Individual
function print_indi(entry::Individual)
  println("The following are the expressions for the individual chromosomes:")
  for i in entry.clist
    println(i.thestring)
  end
  println("This is the header:")
  println(entry.header)
  println("Resulting in an individual with expression:")
  println(entry.thestring)
  println("With repect to the...")
  println("de: ", entry.de, "; bc: ", entry.bc, "; ival: ", entry.ival)
  println("The fitness is:")
  println(entry.fitness)
  println("Divided as:")
  println("error: ", entry.error, "; penalty: ", entry.penalty, "; shape error: ", entry.shape_error)
end

# Julia interprets juxtaposition of a floatingpoint number and a variable as
# multiplication, however the cone here works on explicit multiplication signs
# This is to remove those occurences
function safe_string(instring::String)
  # regex from http://www.regular-expressions.info/floatingpoint.html
  regex = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?x"
  while typeof(match(regex, instring)) != Void
    m = match(regex, instring).match[1:end-1]
    instring = replace(instring, regex, m * " * x", 1)
  end
  return instring
end

# Function to resolve mathematical expressions. <e> symbolises the expression,
# <ep> its first derivative, <epp...> higher derivatives.
function parse_expr(expr::String, indi::Individual)
  name = "e"
  thestring = indi.thestring
  expr = replace(expr, "<" * name * ">", thestring)
  while typeof(match(r"(?<=(\<)).*?(?=\>)", expr)) != Void
    name = name * "p"
    try
      #XXX: You will be able to get an infinite loop here unless you use
      # Calculus with https://github.com/johnmyleswhite/Calculus.jl/pull/107
      # merged
      thestring = differentiate(thestring, :x)
      thestring = string(thestring)
    catch
      thestring = "Inf"
    end
    expr = replace(expr, "<" * name * ">", thestring)
  end
  return expr
end

# Function to initialise a full, e.g. diff-eq dep., Individual
function init_full_indi(inindi::Individual, de::String, bc::Array{Any,1}, ival::Array{Float64,1})
  indi = deepcopy(inindi)
  # Calculate error
  eval_de = parse_expr(de, indi)
  xlist = linspace(ival[1],ival[2],10)
  error = 0
  for x in xlist
    try
      #FIXME: Find a better work-around for the below
      work = eval_de
      work = replace(work, "exp", "å")
      work = replace(work, "x", "$x")
      work = replace(work, "å", "exp")
      work = safe_string(work)
      error = error + (eval(parse(work)))^2
    catch
      error = Inf
    end
  end

  # Calculate penalty
  eval_bc = parse_expr(bc[1], indi)
  penalty = Inf
  try
    work = eval_bc
    work = replace(work, "exp", "å")
    work = replace(work, "x", "$bc[2]")
    work = replace(work, "å", "exp")
    work = safe_string(work)
    penalty = (eval(parse(work)))^2
  catch
    penalty = Inf
  end

  # Calculate shape error
  shape_error = 0
  #= XXX: Disabled
  # This is disabled for now. Theory here is to give a more accurate
  # fitness-value, e.g. a function differing by only a constant is quite close,
  # by a linear term, also quite good, etc, however it is not currectly designed
  # correctly. Was thinking \sum_p \lambda^p \partial_x^p (de), with \lambda
  # a 1/decay factor, but it not right.
  depth = 1   #TODO: Loop the below, and let level determine how far the shape should go
  decay = 100000.0
  try
    eval_de = string(differentiate(eval_de, :x))
  catch
    eval_de = Inf
  end
  shape_error = 0
  for x in xlist
    try
      #FIXME: Find a better work-around for the below
      work = eval_de
      work = replace(work, "exp", "å")
      work = replace(work, "x", "$x")
      work = replace(work, "å", "exp")
      work = safe_string(work)
      shape_error = shape_error + 1/decay*(eval(parse(work)))^2
    catch
      shape_error = Inf
    end
  end
  shape_error = 0
  =#

  # fitness
  pe_bal = 1
  fitness = error + pe_bal*penalty + shape_error

  return Individual(indi.clist, indi.header, indi.thestring, fitness, error, penalty, shape_error, de, bc, ival, indi.head, indi.head_l, indi.tail, indi.tail_l, indi.dict)
end

# When a Chromosome has gotten a new clist, the rest of the entries
# should be reevaluated, which this function takes care of
function reparse_chromo(chromo::Chromosome)
  thestring = parse_clist(chromo.clist, chromo.dict)
  tree = parse_tree(chromo.clist)
  return Chromosome(chromo.clist, thestring, tree, chromo.head, chromo.head_l, chromo.tail, chromo.tail_l, chromo.dict)
end

# When a clist of a Chromosome has changed, this need to act on the individual
function reparse_indi(inindi::Individual)
  indi = deepcopy(inindi)
  new_clist = Chromosome[]
  change = false
  for mem in indi.clist
    if mem.clist != replace(replace(mem.tree,"(",""),")","")  #XXX: This won't work. Tree is truncated compared to clist!
      push!(new_clist, reparse_chromo(mem))
      change = true
    else
      push!(new_clist, mem)
    end
  end

  # Make sure that init_full_indi() only runs when it has indeed changed (see XXX above)
  if change
    indi.clist = new_clist
    thestring = ""

    for i in (length(new_clist)-1):-1:1
      if string(indi.header[i]) != "z"
        thestring = string(indi.header[i]) * "(" * indi.clist[i+1].thestring * ")" * thestring
      end
    end
    thestring = "(" * indi.clist[1].thestring * ")" * thestring

    indi.thestring = thestring
    return init_full_indi(indi, indi.de, indi.bc, indi.ival)
  else
    return indi
  end
end

# Determines if two Individuals are equal
function eq_indi(indi1, indi2)
  eq = true
  for i in 1:length(indi1.clist)
    if !(indi1.clist[i].tree == indi2.clist[i].tree)
      eq = eq && false
    end
  end
  return eq
end

# Generates a population
function gen_pop(pop_size::Int, de::String, bc::Array{Any,1}, ival::Array{Float64,1}, len::Int=1, head::Array{String,1}=head, head_l::Int=head_l, tail::Array{String,1}=tail, tail_l::Int=tail_l, dict::Dict=dict)
  pop = Individual[]
  while length(pop) < pop_size
    indi = init_flat_indi(len, head, head_l, tail, tail_l, dict)
    indi = init_full_indi(indi, de, bc, ival)
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

# Sort population w.r.t. fitness.
function sort_pop(pop::Array{Individual,1})
  return sort(pop, by=x -> x.fitness)
end

### GP Operators ###

# Mutate -- Chromosome level

## Deps for mutate

# Takes an Chromosome and returns a new one
# Changes entry to one of the same type
function mut_change(chromo::Chromosome, mrate::Float64)
  new_clist = ""
  for part in chromo.clist
    if rand() <= mrate
      if string(part) in functions
        new_clist *= rand(functions)
      elseif string(part) in operators
        new_clist *= rand(operators)
      elseif string(part) in terminators
        new_clist *= rand(terminators)
      else
        println("WARNING: '$part' not found in the lists 'functions', 'operators', nor 'terminators'. Nothing done.")
        new_clist *= string(part)
      end
    else
      new_clist *= string(part)
    end
  end
  # thestring, and tree, are left for reparse_indi() to take care of
  return Chromosome(new_clist, chromo.thestring, chromo.tree, chromo.head, chromo.head_l, chromo.tail, chromo.tail_l, chromo.dict)
end

# Takes an Chromosome and returns a new one
# Randomizes entries
function mut_random(chromo::Chromosome, mrate::Float64)
  new_clist = ""
  iter = 1
  for part in chromo.clist
    if rand() <= mrate
      if iter <= chromo.head_l
        new_clist *= rand(chromo.head)
      else
        new_clist *= rand(chromo.tail)
      end
    else
      new_clist *= string(part)
    end
    iter += 1
  end
  # thestring, and tree, are left for reparse_indi() to take care of
  return Chromosome(new_clist, chromo.thestring, chromo.tree, chromo.head, chromo.head_l, chromo.tail, chromo.tail_l, chromo.dict)
end

# Truncates branches
function mut_trunc(chromoin::Chromosome)
  chromo = deepcopy(chromoin)

  # Truncation only acts on chromos with more than one branch
  if length(chromo.tree) == 3
    return chromo
  end

  branches = find( x -> x == '(', chromo.tree)

  # Select a random branch
  b_start = rand(branches)
  # construct the branch
  i = b_start+1
  level = 1
  branch = "("
  while level != 0
    entry = chromo.tree[i]
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
    return chromo
  end

  # flatten it
  branch = replace(replace(branch,"(",""),")","")

  #XXX: Note that this might cause unwanted replacements
  new_clist = replace(chromo.clist, branch, rand(chromo.tail))

  # Pad it to cannonical length
  while length(new_clist) < length(chromo.clist)
    new_clist *= string(rand(chromo.tail))
  end

  return Chromosome(new_clist, chromo.thestring, chromo.tree, chromo.head, chromo.head_l, chromo.tail, chromo.tail_l, chromo.dict)
end

# Grows branches out of terminators
function mut_grow(chromoin::Chromosome)
  chromo = deepcopy(chromoin)

  # Find a terminator that is to be replaced by a tree
  termins = find( x -> string(x) in chromo.tail, chromo.tree)
  pos = rand(termins)

  # Make a new tree
  #XXX: Can be a long loop... (but seems to work ok)
  new_tree = ""
  while length(new_tree) <= 3
    new_tree = init_clist(chromo.head, chromo.head_l, chromo.tail, chromo.tail_l)
    new_tree = parse_tree(new_tree)
  end

  # Insert new tree
  new_clist = String[]
  for i in 1:length(chromo.tree)
    if i == pos
      push!(new_clist, new_tree)
    else
      push!(new_clist, string(chromo.tree[i]))
    end
  end
  new_clist = string(new_clist...)
  new_clist = replace(replace(new_clist,"(",""),")","")

  # complete the clist
  old_clist = replace(replace(chromo.tree,"(",""),")","")
  new_clist = replace(chromo.clist, old_clist, new_clist)

  # Make cannoncical length
  new_clist = new_clist[1:length(chromo.clist)]

  # Make it safe
  nnew_clist = ""
  for i in 1:length(new_clist)
    char = string(new_clist[i])
    if char in chromo.head && !(char in chromo.tail) && i > chromo.head_l
      nnew_clist *= rand(chromo.tail)
    else
      nnew_clist *= string(new_clist[i])
    end
  end
  new_clist = nnew_clist

  # return
  return Chromosome(new_clist, chromo.thestring, chromo.tree, chromo.head, chromo.head_l, chromo.tail, chromo.tail_l, chromo.dict)
end

#TODO
# Swaps arguments for operators
function mut_swap(chromo::Chromosome)
  # Identify operator
  # identify the two following sub-branches
  # swap them
  return chromo
end

## Actual mutate
#TODO: Finish the other methods
#TODO: Some meta-programming in the below should reduce the code-length significantly e.g. eval(parse("mut_" * method))
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
    return indi
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

# Mutate -- Individual level

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
    return p1, p2
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
