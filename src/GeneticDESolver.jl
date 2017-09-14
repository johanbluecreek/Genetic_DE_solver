using Calculus
#XXX: The above cause warnings if running this file with @everywhere, but seems to work...

using Iterators

import Base.show

# Types
include("types.jl")
# General functions
include("genfunc.jl")
# Genetic functions and operators
## population
include("genops/population.jl")
## mutation
include("genops/mutation.jl")
## mutation of head
include("genops/muthead.jl")
## crossover
include("genops/crossover.jl")
## parent selection
include("genops/parentsel.jl")

###########################################
# XXX OLD STUFF XXX
###########################################

### GP Operators ###

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
