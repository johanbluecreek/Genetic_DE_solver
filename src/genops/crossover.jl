############
# Crossover functions
############

# Description
# Main function is
#  crossover: acts on Chromosomes or Individuals, returns (copy of) the same.
# crossover calls specific method functions
#  cross_$method
#

"""
    cross_safe1point(chromo1, chromo2, gselchance)

Performs a safe one-point crossover between two chromosomes.

`gselchance` determines wheter a Gene in the gene-lists of the Chromosomes should be
selected for the crossover operation. When selected, a one-point crossover will take place
that preserves the tree-structure ("safe").

It should be noted that this means that Genes with only one active entry, that is, starting
with a terminator, will not be affected. Genes with more active entries can be affected
(depending on `gselchance`) after the first entry.

# Examples

## Simple example
```julia-repl
julia> include("./src/GeneticDESolver.jl"); include("./src/testenv.jl");

julia> test_chromo1, test_chromo2
((x)/(log((2)/(6)))/(pi)/(7)*(1),(sin(1))+(9)+(6)-(x))

julia> cross_safe1point(test_chromo1, test_chromo2)
((x)/(log((0)-(((pi)-((x)-(2)))/((3)+(5)))))/(pi)/(7)*(1),(sin(0))+(9)+(6)-(x))
```

In the above example we can see that the first Gene was selected: The first Gene in
`test_chromo2` changed from `sin(1)` to `sin(0)`. On the level of entry-lists, where the
function acts, the one-point crossover changed the first Gene's `elist` according to
```
|              | before                         | after                          |
|--------------|--------------------------------|--------------------------------|
| test_chromo1 | x0ucy*98++788p9x09021098947000 | x1p5y02yx7x810562935py5x5x1y11 |
| test_chromo2 | s1p5y02yx7x810562935py5x5x1y11 | s0ucy*98++788p9x09021098947000 |
```

## Advanced example

If we let the Chromosomes grow a bit such that there are more active entries, the crossover
becomes less trivial
 ```julia-repl
julia> test_chromo3 = mutate(mutate(test_chromo1, 1.0, 1.0, "grow"), 1.0, 1.0, "grow")
((log((sin(y))))*(y))/(log(((2)+(cos(x)))/(log(y))))/((8)*(sin((cos(sin(0)))*(cos(7)))))/(sin(log(4)))*(sin((8)/(cos(((log(pi))-(y))-(9)))))

julia> test_chromo4 = mutate(mutate(test_chromo2, 1.0, 1.0, "grow"), 1.0, 1.0, "grow")
(sin(cos((pi))))+((cos(x))/(((6)*(y))/(cos(exp(log(9))))))+((2)*(sin(1)))-((x)+(((log(x))*(x))+(((x)+(5))+(3))))

julia> cross_safe1point(test_chromo3, test_chromo4)
(((log((pi)))*(pi))/(log(((2)+((6)*(y)))/(cos(exp(log(9))))))/((8)*(sin(1)))/(sin(x))*(sin((8)/(x))),
(sin(cos((sin(y)))))+((cos(x))/((cos(x))/(log(y))))+((2)*(sin((cos(sin(0)))*(cos(7)))))-((log(4))+(1)))
```

Where the same table as in the previous example (that is for the first Gene) looks like
```
|              | before                         | after                          |
|--------------|--------------------------------|--------------------------------|
| test_chromo1 | *lusyy0ucyx987p788p947y0902109 | *lupp5y02yx7x810562935py5x5x1y |
| test_chromo2 | scupp5y02yx7x810562935py5x5x1y | scusyy0ucyx987p788p947y0902109 |
```
and we see that the `(pi)` branch has been substituted for `(sin(y))`.
"""
function cross_safe1point(chromoin1::Chromosome, chromoin2::Chromosome, gselchance::Float64=0.8)
    chromo1, chromo2 = deepcopy(chromoin1), deepcopy(chromoin2)
    for g in 1:length(chromo1.glist)
        if rand() <= gselchance
            l1 = length(replace(replace(chromo1.glist[g].tree,"(",""),")",""))
            l2 = length(replace(replace(chromo2.glist[g].tree,"(",""),")",""))
            rnd = 1
            if l1 <= l2
                rnd = rand(1:l1)
            elseif l2 < l1
                rnd = rand(1:l2)
            else
                rnd = rand(2:Int(length(chromo1.glist[1].elist)))
            end
            new_elist1 = chromo1.glist[g].elist[1:rnd] * chromo2.glist[g].elist[rnd+1:end]
            new_elist2 = chromo2.glist[g].elist[1:rnd] * chromo1.glist[g].elist[rnd+1:end]
            chromo1.glist[g].elist = new_elist1
            chromo2.glist[g].elist = new_elist2
            chromo1.glist[g] = reparse_gene(chromo1.glist[g])
            chromo2.glist[g] = reparse_gene(chromo2.glist[g])
        end
    end
    return reparse_chromo(chromo1), reparse_chromo(chromo2)
end

function cross_safe2point(chromoin1::Chromosome, chromoin2::Chromosome, gselchance::Float64=0.8)
    #TODO: Implement this
    chromo1, chromo2 = deepcopy(chromoin1), deepcopy(chromoin2)
    println("WARN: `cross_safe2point()` not implemented yet. Returning input.")
    return chromo1, chromo2
end

"""
    cross_1point(chromo1, chromo2, gselchance)

Unsafe/naïve one-point crossover between two Chromosomes.

The unsafe one-point crossover makes a crossover at any point of the `elist` of a `Gene`.
This should be compared with `cross_safe1point()` which only pickes points for crossover in
active `elist` entries.

This is "unsafe" since the Chromosomes returned are likely to only have exchanged non-active
genetic matter which may cause problems for conserving diversity in a population.
"""
function cross_1point(chromoin1::Chromosome, chromoin2::Chromosome, gselchance::Float64=0.8)
    chromo1, chromo2 = deepcopy(chromoin1), deepcopy(chromoin2)
    for g in 1:length(chromo1.glist)
        if rand() <= gselchance
            rnd = rand(1:length(chromo1.glist[1].elist))
            new_elist1 = (chromo1.glist[g]).elist[1:rnd] * (chromo2.glist[g]).elist[(rnd+1):end]
            new_elist2 = (chromo2.glist[g]).elist[1:rnd] * (chromo1.glist[g]).elist[(rnd+1):end]
            chromo1.glist[g].elist = new_elist1
            chromo2.glist[g].elist = new_elist2
            chromo1.glist[g] = reparse_gene(chromo1.glist[g])
            chromo2.glist[g] = reparse_gene(chromo2.glist[g])
        end
    end
    return reparse_chromo(chromo1), reparse_chromo(chromo2)
end

"""
    cross_2point(chromo1, chromo2, gselchance)

Unsafe/naïve two-point crossover between two Chromosomes.
"""
function cross_2point(chromoin1::Chromosome, chromoin2::Chromosome, gselchance::Float64=0.8)
    chromo1, chromo2 = deepcopy(chromoin1), deepcopy(chromoin2)
    for g in 1:length(chromo1.glist)
        if rand() <= gselchance
            rnd1 = rand(1:length(chromo1.glist[1].elist))
            rnd2 = rand(1:length(chromo1.glist[1].elist))
            while rnd1 == rnd2
                rnd2 = rand(1:length(chromo1.glist[1].elist))
            end
            if rnd1 > rnd2
                rnd1, rnd2 = rnd2, rnd1
            end
            new_elist1 = chromo1.glist[g].elist[1:rnd1] * chromo2.glist[g].elist[(rnd1+1):rnd2] * chromo1.glist[g].elist[(rnd2+1):end]
            new_elist2 = chromo2.glist[g].elist[1:rnd1] * chromo1.glist[g].elist[(rnd1+1):rnd2] * chromo2.glist[g].elist[(rnd2+1):end]
            chromo1.glist[g].elist = new_elist1
            chromo2.glist[g].elist = new_elist2
            chromo1.glist[g] = reparse_gene(chromo1.glist[g])
            chromo2.glist[g] = reparse_gene(chromo2.glist[g])
        end
    end
    return reparse_chromo(chromo1), reparse_chromo(chromo2)
end

"""
    cross_unsafe1point(chromo1, chromo2, gselchance)

Selects one active point of the two chomosomes and mess them up.

This is here to emulate a certain property of the BNF-grammar. In the BNF, where this code's
version of the `elist` of a chromosome would be a string of numbers, what those numbers
represent; operator, function or digit, depends on all previous entries. This differs from
the Polish-notation used here, an `s` in the elist is always a `sin`-function (as per the
default dictionary used here). This means in a one-point crossover in BNF, while still
preserving the genetic material, the genetic material will in general get a new mathematical
expressions.

We cannot have both here. If we preserve genetic material the mathematical expression will
be the same. This function emulates the occurence of random behaviour of the genetic
material used at the cost of preserving genetic material.

Similar to this there is also the `unsaferandom` method for `mutate()`.
"""
function cross_unsafe1point(chromoin1::Chromosome, chromoin2::Chromosome, gselchance::Float64=0.8)
    chromo1, chromo2 = deepcopy(chromoin1), deepcopy(chromoin2)
    for g in 1:length(chromo1.glist)
        if rand() <= gselchance
            l1 = length(replace(replace(chromo1.glist[g].tree,"(",""),")",""))
            l2 = length(replace(replace(chromo2.glist[g].tree,"(",""),")",""))
            rnd = 1
            if l1 <= l2
                rnd = rand(1:l1)
            elseif l2 < l1
                rnd = rand(1:l2)
            else
                rnd = rand(2:Int(length(chromo1.glist[1].elist)))
            end
            #XXX: init_elist() uses global variables, that are not passed to this function
            # so that they could be overwritten! Rethink global variables!
            # Either use them in function arguments, but then always so! or never! Or use
            # the types for this!
            etail1 = init_elist()[rnd+1:end]
            etail2 = init_elist()[rnd+1:end]
            new_elist1 = chromo1.glist[g].elist[1:rnd] * etail1
            new_elist2 = chromo2.glist[g].elist[1:rnd] * etail2
            chromo1.glist[g].elist = new_elist1
            chromo2.glist[g].elist = new_elist2
            chromo1.glist[g] = reparse_gene(chromo1.glist[g])
            chromo2.glist[g] = reparse_gene(chromo2.glist[g])
        end
    end
    return reparse_chromo(chromo1), reparse_chromo(chromo2)
end

"""
    cross_random(chromo1, chromo2, gselchance)

Random jumps of Genes between two Chromosomes.

Instead of acting on the `elist` of Genes, the whole Genes are transfered between the
Chromosomes. `gselchance` determines a chance for the Genes of `chromo2` to replace one of
Genes in `chromo1`. That is `gselchance = 0.5` will randomize the genes, `gselchance = 0.1`
will preserve most genes of `chromo1`.
"""
function cross_random(chromoin1::Chromosome, chromoin2::Chromosome, gselchance::Float64=0.5)
    chromo1, chromo2 = deepcopy(chromoin1), deepcopy(chromoin2)
    new_glist1 = Gene[]
    new_glist2 = Gene[]
    for g in 1:length(chromo1.glist)
        if rand() < gselchance
            push!(new_glist1, chromo1.glist[g])
            push!(new_glist2, chromo2.glist[g])
        else
            push!(new_glist1, chromo2.glist[g])
            push!(new_glist2, chromo1.glist[g])
        end
    end
    chromo1.glist = new_glist1
    chromo2.glist = new_glist2
    return reparse_chromo(chromo1), reparse_chromo(chromo2)
end


"""
    crossover(chromo1, chromo2, gselchance, method)

Performs crossover between two Chromosomes using `method`.

# Methods

Avaliable methods are

`["safe1point", "1point", "2point", "random"]`

See documentation for `cross_[method]` for a description for each of the methods.

# Examples

## Random
```julia-repl
julia> include("./src/GeneticDESolver.jl"); include("./src/testenv.jl");

julia> test_chromo1, test_chromo2
((y)*((2)+(exp(log(1))))*((0)+(y))*(x),(6)/(4)*(0)+(y)/(sin(8)))

julia> crossover(test_chromo1, test_chromo2, 0.5, "random")
((y)*((2)+(exp(log(1))))*(0)*(x),(6)/(4)*((0)+(y))+(y)/(sin(8)))
```

## Two-point

```julia-repl
julia> include("./src/GeneticDESolver.jl"); include("./src/testenv.jl");

julia> test_chromo1, test_chromo2
((y)*(x)-(x)-(x)+((x)-(exp(y))),((pi)+(x))+(sin(7))*(cos(cos(x)))*(2)*(y))

julia> map(x -> x.elist, test_chromo1.glist), map(x -> x.elist, test_chromo2.glist)
(
    String["yxsxuy15xx20y0xp6y9y9pp0x40631","xcxl*9y6u254761y965x9169872870","x3pyxy86p7x32y3y38x6y71p7xyx0p","x1x-y2l/3u198yp36p60x99926446y","-xey31yu4y14414p21yyx69p80yx43"],
    String["+pxux8yc7u27361x11y42x68447530","s7+5sx/x*y0y15249798pp89p93270","ccxyl-uye-771p1p4028562574x1x7","2x**8*lxyxx0y7y5y73y20x52yp14x","y00cx4xl9x775y4p09862x8015p8yp"]
)

julia> c3, c4 = crossover(test_chromo1, test_chromo2, 0.8, "2point")
((y)*(x)-(x)-(x)+((x)-(exp(y))),((pi)+(x))+(sin(7))*(cos(cos(x)))*(2)*(y))

julia> map(x -> x.elist, c3.glist), map(x -> x.elist, c4.glist)
(
    String["yxsxuy1c7u27361x11y42x68447531","xcxl*9y6u254761y965x9169872870","x3pyl-uye-771p1p4028571p7xyx0p","x1x-8*lxyxx0y7y5yp60x99926446y","-xey31yl9x775y4p09862x8015p843"],
    String["+pxux8y5xx20y0xp6y9y9pp0x40630","s7+5sx/x*y0y15249798pp89p93270","ccxyxy86p7x32y3y38x6y62574x1x7","2x**y2l/3u198yp3673y20x52yp14x","y00cx4xu4y14414p21yyx69p80yxyp"]
)
```

"""
function crossover(chromoin1::Chromosome, chromoin2::Chromosome, gselchance::Float64=0.8, method = "safe1point")
    chromo1, chromo2 = deepcopy(chromoin1), deepcopy(chromoin2)
    methods = ["safe1point", "1point", "2point", "unsafe1point", "random"]
    if method in methods
        ee = "cross_$method"
        ee = eval(parse(ee))
        return ee(chromo1, chromo2, gselchance)
    else
        println("WARNING: No support for method: '$method'. Nothing done.")
        return chromo1, chromo2
    end
end

function crossover(indiin1::Individual, indiin2::Individual, gselchance::Float64=0.8, method = "safe1point", cselchance::Float64=0.8)
    indi1, indi2 = deepcopy(indiin1), deepcopy(indiin2)
    new_clist1 = Chromosome[]
    new_clist2 = Chromosome[]
    for c in 1:length(indi1.clist)
        if rand() < cselchance
            c1, c2 = crossover(indi1.clist[c], indi2.clist[c], gselchance, method)
            push!(new_clist1, c1)
            push!(new_clist2, c2)
        else
            push!(new_clist1, indi1.clist[c])
            push!(new_clist2, indi2.clist[c])
        end
    end
    indi1.clist = new_clist1
    indi2.clist = new_clist2
    return reparse_indi(indi1), reparse_indi(indi2)
end
