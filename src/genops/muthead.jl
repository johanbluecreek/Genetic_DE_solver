############
# Mutation of heads
############

# Description
# Main function is
#  muthead: acts on Chromosomes or Individuals, returns (copy of) the same.
# muthead calls specific method functions
#  muth_$method: acts on Chromosomes, returns (copy of) the same.
#

"""
    muth_scramble(chromo::Chromosome[, mrate])

Mutates the head of `chromo` by scrambling the order of Genes.

`muth_scramble()` has an optional argument `mrate` that has no effect, and is present only
for uniform `muth_method()` syntax.

# Examples
```julia-repl
julia> include("./src/GeneticDESolver.jl");

julia> include("src/testenv.jl");

julia> test_chromo
(log((6)*(x)))+(y)/(2)+(x)/(x)

julia> muth_scramble(test_chromo)
(x)+(2)/(x)+(log((6)*(x)))/(y)
```
"""
function muth_scramble(inchromo::Chromosome, mrate::Float64=0.6)
    chromo = deepcopy(inchromo)
    if length(chromo.glist) > 1
        #XXX: Does this create copies of the Genes or not?
        new_glist = map(x->chromo.glist[x], randperm(length(chromo.glist)))
        chromo.glist = new_glist
    else
        #XXX: Make proper warning calls of these
        println("WARNING: Not enough genes to perform 'scramble'")
    end
    return reparse_chromo(chromo)
end

"""
    muth_jump(chromo::Chromosome[, mrate])

Mutates the head of `chromo` by making two Genes jump places.

`muth_jump()` has an optional argument `mrate` that has no effect, and is present only
for uniform `muth_method()` syntax.

# Examples
```julia-repl
julia> include("./src/GeneticDESolver.jl");

julia> include("src/testenv.jl");

julia> test_chromo
(log((6)*(x)))+(y)/(2)+(x)/(x)

julia> muth_jump(test_chromo)
(2)+(y)/(log((6)*(x)))+(x)/(x)
```
"""
function muth_jump(inchromo::Chromosome, mrate::Float64=0.6)
    chromo = deepcopy(inchromo)
    if length(chromo.glist) > 1
        r1 = rand(1:length(chromo.glist))
        r2 = rand(1:length(chromo.glist))
        while r1 == r2
            r2 = rand(1:length(chromo.glist))
        end
        if r1 > r2
            r1, r2 = r2, r1
        end
        new_glist = Gene[]
        append!(new_glist, chromo.glist[1:(r1-1)])
        append!(new_glist, chromo.glist[r2:r2])
        append!(new_glist, chromo.glist[r1+1:(r2-1)])
        append!(new_glist, chromo.glist[r1:r1])
        append!(new_glist, chromo.glist[r2+1:end])
        chromo.glist = new_glist
    else
        println("WARNING: Not enough genes to perform 'jump'")
    end
    return reparse_chromo(chromo)
end

"""
    muth_combo(chromo::Chromosome[, mrate])

Mutates the head of `chromo` by changing how the genes are combined.

`muth_combo()` has an optional argument `mrate` which sets the chance a head operator in
the Chromosome will change randomly.

# Examples
```julia-repl
julia> include("./src/GeneticDESolver.jl");

julia> include("src/testenv.jl");

julia> test_chromo
(log((6)*(x)))+(y)/(2)+(x)/(x)

julia> muth_combo(test_chromo, .6)
(log((6)*(x)))*(y)*(2)+(x)*(x)
```
"""
function muth_combo(inchromo::Chromosome, mrate::Float64=0.6)
    chromo = deepcopy(inchromo)
    new_header = ""
    for part in chromo.header
      if rand() <= mrate
        new_header *= rand(operators)
      else
        new_header *= string(part)
      end
    end
    chromo.header = new_header
    return reparse_chromo(chromo)
end

"""
    muthead(input, mrate, method)

Mutates `input` on the level of the head (that is, how `Gene`s are combined).

# Methods

`["scramble", "jump", "combo"]`

See documentation for `muth_[method]` for a description for each of the methods.

# Examples
```julia-repl
julia> include("./src/GeneticDESolver.jl");

julia> include("src/testenv.jl");

julia> test_chromo
(log((6)*(x)))+(y)/(2)+(x)/(x)

julia> muthead(test_chromo, .6, "jump")
(2)+(y)/(log((6)*(x)))+(x)/(x)

julia> muthead(test_chromo, .6, "scramble")
(x)+(2)/(x)+(log((6)*(x)))/(y)

julia> muthead(test_chromo, .6, "combo")
(log((6)*(x)))*(y)*(2)+(x)*(x)
```
"""
function muthead(inchromo::Chromosome, mrate::Float64=0.6, method::String="jump")
    chromo = deepcopy(inchromo)
    methods = ["scramble", "jump", "combo"]
    if method in methods
        ee = "muth_$method"
        ee = eval(parse(ee))
        ee = ee(chromo, mrate)
        return ee
    else
        println("WARNING: No support for method: '$method'. Nothing done.")
        return chromo
    end
end
