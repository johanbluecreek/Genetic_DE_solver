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
