
using ProgressMeter

# Include all functions
include("../src/GeneticDESolver.jl")
# this is defined outside of `TL_Solver` to not screw up world age.

function TL_solver(de_file::String, penalty_factor::Union{Int64,Float64}=100, shape_decay::Union{Int64,Float64}=Inf, pop_size::Int64=1000, stop::Int64=2000, sens::Float64=10.0^(-7))

    # Include the differential equation
    include(de_file)

    # Redefine the global defaults
    global FLIST = flist
    global VARS = vars

    # Set function and operator basis
    functions = ["s", "c", "e", "l"];
    operators = ["+", "-", "*", "/"];
    # as well as digits
    digits = vcat(["$i" for i=range(0,10)]);
    # These are as in [TL].

    # Redefine the global defaults
    global FUNCTIONS = functions
    global OPERATORS = operators
    global DIGITS = digits

    # Let the functions know which are the
    terminators = vcat(digits, vars);

    # Redefine the global defaults
    global TERMINATORS = terminators

    # The grammar of [TL] is different,
    # Chance to select any operator is equal to the chance of getting the unity operator,
    # any function, or any digit. However, the variables x, y, z also have that same chance.
    # Meaning vcat(digit, vars) will not be correct since this would give, e.g., '0' the same
    # chance of being selected as 'x' -- which is not what [TL] has. The correct balance between
    # digits and variables would be vcat(digit, repeat(vars, outer=10)) making the variable
    # as likely as any digit.
    # This means, normalising using vars, or the unity operator, the correct head is
    head = vcat(
        repeat(operators, outer=5),
        repeat(["u"], outer=20),
        repeat(functions, outer=5),
        repeat(digits, outer=2),
        repeat(vars, outer=20));
    # This is not completely correct, since the random number chosen by [TL] is 0:255, which
    # will give an uneven distribution (due to modulus operations). I will not try and correct
    # for this here.
    tail = vcat(digits, repeat(vars, outer=10));
    # [TL] states that they use fixed length, although they do not specify how they fix the
    # length, here we take the fixed length to be in terms of head_l, and we set the tail to be
    # made out of head, having removed operators, unitary operator, and functions.
    head_l = 50;
    tail_l = 51;

    # Redefine the global defaults
    global HEAD = head
    global TAIL = tail
    global HEAD_L = head_l
    global TAIL_L = tail_l

    # [TL] does not have their expression composed of several genes in this sense. But they do
    # have "wrapping events" that seems to be occuring when maximum length is reached. Hence
    # we will simulate this with allowing for two genes, but mostly 'z' operators allowed for
    # the header (which kills one gene, hence no wrapping)
    glen = 2
    header_operators = vcat(operators, repeat(["z"], outer=20));
    #XXX: Note that this is not too good. If length of the active elements of a gene were
    # visible, then a mutation of the header could be triggered by that.

    # Redefine global defaults
    global GLEN = glen
    global HEADER_OPERATORS = header_operators

    # [TL] has a penalty factor for the boundary conditions contribution to the complete
    # fitness. They also do not have a shape implemented, so we disable that here. These are
    # given as an argument to the function, so we simply overide defaults
    global PENALTY_FACTOR = penalty_factor
    global SHAPE_DECAY = shape_decay

    pop = gen_pop(pop_size, de, bc, ival)

    #Allow only sane expressions:
    for i in 1:length(pop)
        if (isequal(pop[i].fitness, NaN) || isequal(pop[i].fitness, Inf))
            while (isequal(pop[i].fitness, NaN) || isequal(pop[i].fitness, Inf))
                pop[i] = gen_pop(1, de, bc, ival)[1]
            end
        end
    end

    pop = sort_pop(pop)

    iter = 1
    thetime = time()
    progIter = Progress(stop, 0.1, de_file)


    # Start solving
    while iter < stop && pop[1].fitness > sens

        # [TL] uses only two genetic operations:
        #  crossover: one-point (here it is taken to be the 'safe' one-point crossover)
        #  mutation: random

        # Parents are selected by tournament:

        # Number of pair of parents:
        p_pop_size = Int(round((1-0.1)*pop_size/2))
        # [TL] calls (twice) this c = (1-s)*g, with: g = pop_size, s = replication rate (set by
        # [TL] s = 0.1) hence c = (1-.1)*1000 = 900.
        p_pop = repeat(gen_p_pop(pop, 1, rand(2:pop_size), "tournament"), outer=p_pop_size)
        # [TL] uses 'tournamnet' of a 'K' sized randomly selected sub-pop. [TL] does not specify
        # K more than 'K ≥ 2'. Assuming K is randomly chosen here, and not fixed, the above is
        # split to have one new K for each parent pair.

        # Parents selected, compute the crossover:
        half = Int(round(p_pop_size))
        c_pop = map(x -> crossover(x..., 1.0, "unsafe1point", 1.0), p_pop[1:half])
        c_pop = vcat(c_pop, map(x -> crossover(x..., 1.0, "safe1point", 1.0), p_pop[half+1:end]))
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
        s_pop = shuffle(pop[1:(pop_size-length(c_pop))])
        half = Int(round(length(s_pop)))
        m_pop = map(x -> mutate(x, 0.05, 1.0, "unsaferandom", 1.0), s_pop[1:half])
        m_pop = vcat(m_pop, map(x -> mutate(x, 0.05, 1.0, "random", 1.0), s_pop[half+1:end]))
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
        # Mutate head to simulate "wrapping"
        pop = map(x -> muthead(x, 0.6, "combo", 0.8), pop)

        # And sort again
        pop = sort_pop(pop)

        iter += 1
        next!(progIter)

    end

    thetime = time()-thetime
    sol = pop[1]
    # Did we actually solve the de? if not...
    if sol.fitness > sens
        sol = Inf
    end

    return [sol, thetime, iter]
end
