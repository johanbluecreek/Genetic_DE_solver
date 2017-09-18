using Calculus
#XXX: The above cause warnings if running this file with @everywhere, but seems to work...

using Iterators

import Base.show

GDES_loaded = true

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

##
# TODO/Comments
#
# * Be careful with the deepcopy() calls! There are alot of unnecessary calls that slows
#   stuff down. (e.g. one in mutate when mut_change has one anyway)
##

###########################################
# XXX OLD STUFF XXX
###########################################

### GP Operators ###


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
