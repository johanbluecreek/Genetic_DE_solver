############
# Mutation functions
############

# Description
# Main function is
#  mutate: acts on Chromosomes or Individuals, returns (copy of) the same.
#   if Chomosome, Gene to be acted on is determined by mselchance
#   if Individual, Chromosomes are selected by cselchance.
# mutate calls specific method functions
#  mut_$method: acts on Genes, returns (copy of) the same.
#

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
function mut_change(genein::Gene, mrate::Float64)
    gene = deepcopy(genein)
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
    gene.elist = new_elist
    return reparse_gene(gene)
end

"""
    mut_random(gene::Gene, mrate::Float64)

Changes entry to a random entry (respecting head/tail).

`mrate` is the likelyhood that an entry will change, e.g. `rand() <= mrate`
determines change.

# Examples
```julia-repl
julia> include("./src/GeneticDESolver.jl");

julia> include("src/testenv.jl");

julia> test_gene
(1)-(((0))+(y))

julia> mut_random(test_gene, .2)
(1)/(((0))+(y))

```
"""
function mut_random(gene::Gene, mrate::Float64=1.0)
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
julia> include("./src/GeneticDESolver.jl");

julia> include("src/testenv.jl");

julia> test_gene
(1)-(((0))+(y))

julia> mut_trunc(test_gene)
(1)-((7)+(y))

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

Since the entry-lists are set to a fixed length (`head_l` + `tail_l`), growing cannot be
done indefinetly in a controled way. Once the growth is taken beyond `head_l`, such that the
expression is no longer garanteed to terminate, `mut_grow()` will remove non-terminators
from the tail part of the entry list, and make it terminate. To avoid running into this
slight erractic behaviour, increase `head_l` and `tail_l` accordingly.

# Examples
```julia-repl
julia> include("./src/GeneticDESolver.jl");

julia> include("src/testenv.jl");

julia> test_gene
(1)-(((0))+(y))

julia> mut_grow(test_gene)
(1)-(((cos((x))))+(y))

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
    nnew_elist = ""
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
    mutate(input, mrate, mselchance, method[, cselchance])

Mutates `input` using `method`.

# Inputs

The `input` can be either `Chromosome` or `Individual`, where for an individual Chromosomes
will be selected for mutation depending on `cselchance`.

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

# See Also

`mutate` mutates `Gene`-expressions. See `muthead` for mutation on the level of
`Chromosomes`.

"""
function mutate(inchromo::Chromosome, mrate::Float64=0.6, mselchance::Float64=0.2, method::String="change")
    chromo = deepcopy(inchromo)
    methods = ["change", "swap", "grow", "trunc", "random"]
    if method in methods
        if length(chromo.glist) == 0 #XXX: Not been tested yet!
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

function mutate(inindi::Individual, mrate::Float64=0.6, mselchance::Float64=0.2, method::String="change", cselchance::Float64=0.8)
    indi = deepcopy(inindi)
    new_clist = Chromosome[]
    for chromo in indi.clist
        if rand() <= cselchance
            push!(new_clist, mutate(chromo, mrate, mselchance, method))
        else
            push!(new_clist, chromo)
        end
    end
    indi.clist = new_clist
    indi = reparse_indi(indi)
    return indi
end
