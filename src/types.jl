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
