############
# Parent selection
############

# Description
# Main function is
#  p_select: Which selects parents for crossover or mutation
# p_select calls specific method functions
#  pmet_$method
#
# See also
#  gen_m_pop
#  gen_p_pop

#XXX: Does the returned populations need to be copies?

"""
    pmet_tournament(pop, selsize)

Parent selection method tournament.

Returns a parent from the population by creating a `selsize` sized sub-population and
choosing the one with the highest fitness.
"""
function pmet_tournament(pop::Array{Individual, 1}, selsize::Int=5)
    pop_size = length(pop)
    p1 = rand(pop,rand(2:selsize))
    p1 = sort_pop(p1)[1]
    return p1
end

"""
    pmet_random(pop, selsize)

Parent selcetion method random.

Returns a parent from the population by selecting a random one. This is simply a short-hand
for `rand(pop)`, and to be used in `p_select`.
"""
function pmet_random(pop::Array{Individual, 1}, selsize::Int=5)
    return rand(pop)
end

"""
    p_select(pop, selsize, method)

Perform a selection of parents from a population `pop` using `method`.

# Methods

`["tournament", "random"]`

See documentation for `pmet_[method]` for a detailed explanation for each method.
"""
function p_select(pop::Array{Individual, 1}, selsize::Int=5, method::String="tournament")
    methods = ["tournament", "random"]
    if method in methods
        ee = "pmet_$method"
        ee = eval(parse(ee))
        return ee(pop, selsize)
    else
        println("WARNING: No support for method: '$method'. Returning first member of population.")
        return pop[1]
    end
end

"""
    gen_p_pop(pop, size, selsize, method)

Generates `size` pair of parents using selection method `method`.

# Examples

```julia-repl
julia> include("./src/GeneticDESolver.jl"); include("./src/testenv.jl");

julia> test_pop = gen_pop(50, test_de, test_bc, test_ival)
50-element Array{Individual,1}:
 2066.999999999999: ((5)-(y)), (y)+((x))
 6027.0: (y)*(2), (x)+((9)*(y))
 1937.7254644660927: (exp(x))+(y), (6)-(x)
 2052.016029981668: ((y)-(cos(x))), (5)*((x))
 200.0: ((x)-(x)), (exp(pi))/(exp(pi))
 ...
 28800.0: (6)*(2), (y)-(y)
 13518.349664073132: (x)/((y)/(6)), (pi)-((pi)*(x))
 14824.484233611513: (log(5)), (2)-((9))
 38777.00000000001: (3)*(((9)-(y))-(5)), ((4)*(x))-(6)

julia> gen_p_pop(test_pop)
10-element Array{Tuple{Individual,Individual},1}:
 (7331.625000000004: (3)/(8), (9)*(y),132.75389417929688: (exp(y))-((x)-(6)), (x)+(7))
 (848.6876340317573: (x)/(y), (sin((5)+(y))),1937.7254644660927: (exp(x))+(y), (6)-(x))
 (101.00000000000001: (x)-(x), (y),11.222222222222214: (0)*(x), (y)/(3))
 (1053.623153167669: (sin(y)), (9)/(pi),3885.0: (x)*(8), (y))
 (0.0: (0)+(x), (x),2205.611833956204: (7)-(x), (pi))
 (195.62367825076157: (x)/(y), (x),411.71178942053984: (x)*(pi), (y))
 (0.0: (0)+(x), (x),5287.195003089709: (1)+((8)-(2)), (pi)+(9))
 (1575.0979087153003: (y), (2)-((1)+(pi)),101.00000000000001: (x)-(x), (y))
 (1885.7246527744217: (x)/(x), (exp(sin(y)))+(2),11.222222222222214: (0)*(x), (y)/(3))
 (0.0: (0)+(x), (x),11.222222222222214: (0)*(x), (y)/(3))
```
"""
function gen_p_pop(pop::Array{Individual, 1}, size::Int=10, selsize::Int=5, method::String="tournament")
    new_pop = Individual[]
    for i in 1:size
        push!(new_pop, p_select(pop, selsize, method))
        next_p = p_select(pop, selsize, method)
        #XXX: This should probably be stopped if it goes on for too long
        j = 1
        while eq_indi(next_p, new_pop[(2*i-1)]) && j <= 5
            next_p = p_select(pop, selsize, method)
            j += 1
        end
        push!(new_pop, next_p)
    end
    return collect(zip(new_pop[1:2:end], new_pop[2:2:end]))
end

"""
    gen_m_pop(pop, size, selsize, method)

Generates `size` number of individuals using selection method `method`.
"""
function gen_m_pop(pop::Array{Individual, 1}, size::Int=10, selsize::Int=5, method::String="tournament")
    new_pop = Individual[]
    for i in 1:size
        push!(new_pop, p_select(pop, selsize, method))
    end
    return new_pop
end
