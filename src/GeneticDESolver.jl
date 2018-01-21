using Calculus
#XXX: The above cause warnings if running this file with @everywhere, but seems to work...

using Iterators

import Base.show

# Variable to determine if this file has been included
GDES_loaded = true

# Define default settings for global variables
include("globaldefault.jl")

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
