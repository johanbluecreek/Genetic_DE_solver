############
# Mutation of heads
############

"""
    muthead(input, mrate, method)

Mutates `input` on the level of the head (that is, how `Gene`s are combined).

# Methods

`["scramble", "jump", "combo"]`

`"scramble"` randomize position of genes; `"jump"` switches places of a pair of genes;
`"combo"` mutates how the genes are combined.

# Examples

"""
function muthead(inchromo::Chromosome, mrate::Float64=0.6, method::String="jump")
    chromo = deepcopy(inchromo)
    methods = ["scramble", "combo"]
    if method in methods
        #XXX: Might want to have separate functions for all of these in the future...
        # See documentation for `muth_[method]` for a description for each of the methods.
        if method == "scramble"
            if length(chromo.glist) > 1
                #XXX: Does this create copies of the Genes or not?
                new_glist = map(x->chromo.glist[x], randperm(length(chromo.glist)))
                chromo.glist = new_glist
            else
                #XXX: Make proper warning calls of these
                println("WARNING: Not enough genes to perform '$method'")
            end
        elseif method == "jump"
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
                append!(new_glist, chromo.glist[r2])
                append!(new_glist, chromo.glist[r1+1:(r2-1)])
                append!(new_glist, chromo.glist[r1])
                append!(new_glist, chromo.glist[r2+1:end])
                chromo.glist = new_glist
            else
                println("WARNING: Not enough genes to perform '$method'")
            end
        elseif method == "combo"
            new_header = ""
            for part in chromo.header
              if rand() <= mrate
                new_header *= rand(operators)
              else
                new_header *= string(part)
              end
            end
            chromo.header = new_header
        end
        chromo = reparse_chromo(chromo)
        return chromo
    else
        println("WARNING: No support for method: '$method'. Nothing done.")
        return chromo
    end
end
