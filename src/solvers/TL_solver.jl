
###
# This is a solver following the specifications of
# [TL] (see README). When executed, it should display
# similar performance to the results of that paper.
###

include("./../GeneticDESolver.jl")
include("./../../odes/ode1.jl")

functions = ["s", "c", "e", "l"];
operators = ["+", "-", "*", "/"];
digits = vcat(["$i" for i=range(0,10)]);

# Even though ode1 only depends on 'x', [TL] uses 'x', 'y' and 'z' in their grammar
# so lets redefine the 'vars' present in 'ode1.jl'
#vars = ["x", "y", "z"]
# which also means we need to redefine ival and bc
#bc = [["(<e>) - 20.1", [("x", 0.1), ("y", Inf), ("z", Inf)]]]
#ival = [(0.1,1.0), (Inf, Inf), (Inf, Inf)]
#XXX: Not using this since it makes it _really_ slow.

# This is just a list of what should be considered terminators. Distribution does not matter
# (that is, vcat(digits, vars) will not be different from vcat(digit, vars, vars))
terminators = vcat(digits, vars);

header_operators = vcat(operators);

# The grammar of [TL] is different,
# Chance to select any operator is equal to the chance of getting the unity operator,
# any function, or any digit. However, the variables x, y, z also have that same chance.
# Meaning vcat(digit, vars) will not be correct since this would give, e.g., '0' the same
# chance of being selected as 'x' -- which is not what [TL] has. The correct balance between
# digits and variables would be vcat(digit, repeat(vars, outer=10)) making the variable
# as likely as any digit.
# This means, normalising using vars, or the unity operator, the correct head is
head = vcat(repeat(operators, outer=5), repeat(["u"], outer=20), repeat(functions, outer=5), repeat(digits, outer=2), repeat(vars, outer=20));
# This is not completely correct, since the random number chosen by [TL] is 0:255, which
# will give an uneven distribution (due to modulus operations). I will not try and correct
# for this here.
tail = vcat(digits, repeat(vars, outer=10));
# [TL] states that they use fixed length, although they do not specify how they fix the
# length, here we take the fixed length to be in terms of head_l, and we set the tail to be
# made out of head, having removed operators, unitary operator, and functions.
head_l = 50;
tail_l = 100;
#XXX: Read about "wrapping events"

# [TL] does not have their expression composed of several genes in this sense.
glen = 1

dict = Dict(
    "+" => "(<expr>)+(<expr>)",
    "-" => "(<expr>)-(<expr>)",
    "*" => "(<expr>)*(<expr>)",
    "/" => "(<expr>)/(<expr>)",
    "s" => "sin(<expr>)",
    "c" => "cos(<expr>)",
    "l" => "log(<expr>)",
    "e" => "exp(<expr>)",
    "u" => "(<expr>)",
);

# [TL] has population size 1000, termination at 2000 itrations or 10.0^(-7) accuracy
pop_size = 1000
stop = 2000
sens = 10.0^(-7)

pop = gen_pop(pop_size, de, bc, ival, flist, glen, header_operators, head, head_l, tail, tail_l, dict)

#Allow only sane expressions:
for i in 1:length(pop)
    if (isequal(pop[i].fitness, NaN) || isequal(pop[i].fitness, Inf))
        while (isequal(pop[i].fitness, NaN) || isequal(pop[i].fitness, Inf))
            pop[i] = gen_pop(1, de, bc, ival, flist, glen, header_operators, head, head_l, tail, tail_l, dict)[1]
        end
    end
end

pop = sort_pop(pop)

iter = 1

start = time()

while iter < stop && pop[1].fitness > sens

    part = time()

    # [TL] uses only two genetic operations:
    #  crossover: one-point (here it is taken to be the 'safe' one-point crossover)
    #  mutation: random

    # Parents are selected by tournament:

    # Number of pair of parents:
    p_pop_size = 450
    # [TL] calls (twice) this c = (1-s)*g, with: g = pop_size, s = replication rate (set by
    # [TL] s = 0.1) hence c = (1-.1)*1000 = 900.
    p_pop = repeat(gen_p_pop(pop, 1, rand(2:pop_size), "tournament"), outer=p_pop_size)
    # [TL] uses 'tournamnet' of a 'K' sized randomly selected sub-pop. [TL] does not specify
    # K more than 'K ≥ 2'. Assuming K is randomly chosen here, and not fixed, the above is
    # split to have one new K for each parent pair.

    # Parents selected, compute the crossover:
    c_pop = map(x -> crossover(x..., 1.0, "safe1point", 1.0), p_pop)
    # g- and cselchance are both selected to be 1.0 since we have only one Gene and one
    # Chromosome per Individual, and it must be selected for crossover (1.0).
    c_pop = vcat(map(x -> x[1], c_pop), map(x -> x[2], c_pop))

    #XXX: Should `NaN` and `Inf` fitness individuals be discarded? If not, set false:
    new_c_pop = Individual[]
    if true
        for i in c_pop
            if !(isequal(i.fitness, NaN) || isequal(i.fitness, Inf))
                push!(new_c_pop, i)
            end
        end
    end
    c_pop = new_c_pop

    # The members of 'c_pop' should replace the worst entries in the population, hence:
    pop = vcat(pop[1:(pop_size-length(c_pop))], c_pop)

    # The mutation is said by [TL] to take place:
    # "The mutation operation is applied to every chromosome excluding those which have been
    # selected for replication in the next generation."
    # which I find slightly ambigous, but I take it to mean the non c_pop members of pop
    m_pop = map(x -> mutate(x, 0.05, 1.0, "random", 1.0), pop[1:(pop_size-length(c_pop))])
    # mrate etc chosen to be similar to [TL].

    #XXX: Should `NaN` and `Inf` fitness individuals be discarded? If not, set false:
    if true
        for i in 1:length(m_pop)
            if (isequal(m_pop[i].fitness, NaN) || isequal(m_pop[i].fitness, Inf))
                while (isequal(m_pop[i].fitness, NaN) || isequal(m_pop[i].fitness, Inf))
                    m_pop[i] = mutate(pop[i], 0.05, 1.0, "random", 1.0)
                end
            end
        end
    end

    # Make the new population
    pop = vcat(m_pop, c_pop)

    # And sort again
    pop = sort_pop(pop)

    println("iter: ", iter)
    println("time: ", time()-part)
    println("best: ", pop[1])
    println("worst: ", pop[end])

    iter += 1


end

println("Done!")
println("Solution: ", pop[1])
println("Found in $iter iterations, and $(time()-start) seconds.")

# Running this you should have something similar as
# Min | Max  | Avg.
#  8  | 1452 | 653
# as the results for the iterations needed.
