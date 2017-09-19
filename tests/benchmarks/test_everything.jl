
# Benchmark search: Search using all available genetic operators

#XXX: WARNING: There seems to be a memory leak that causes this to use >32GB RAM if
# len ~ 20 or so. So try to avoid using that.

###
# Load the package
###
include("./../../src/GeneticDESolver.jl")

###
# Load DE to solve
###
include("./../../odes/ode1.jl")

###
# Set global variables
###

# Basis set
functions = ["s", "c", "e", "l", "u"];
operators = ["+", "-", "*", "/"];
digits = vcat(["$i" for i=range(0,10)], ["p"]);
# Relevant variables for PDEs (can also be used as unspecified constants)

# List of the above that terminates an expression tree
terminators = vcat(digits, vars, vars, vars, vars, vars);

# header_operators combines genes
header_operators = vcat(operators, ["z"]);
# "z" is an operator that deactivates the following chromosome

# Select what should be part of the 'head' of a gene (first head_l entries)
# and which should be part of the 'tail' of a gene (last tail_l entries)
head = vcat(functions, operators, vars, digits, vars, vars, vars, vars);
tail = vcat(digits, vars);

# Leangth of the head and tail of an elist
head_l = 10;
tail_l = 20;

# Number of Genes in a Chromosome
glen = 4

# Dictionary to translate operators and functions to expressions
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
    "p" => "pi"
);

###
# The algorithm
###

# Population size
pop_size = 50
# Number of iterations before stopping
stop = 200
# Sensitivity (numbers smaller than this are "zero")
sens = 10.0^(-7)

# List of times it took to find solution
times = Float64[]
# List of iterations it took to find solution
iters = Int[]
# How many times a solution should be searched for
len = 5


for i = 1:len

    # Iteration counter for the genetic part
    iter = 1
    # Similar for time
    start = time()


    pop = gen_pop(pop_size, de, bc, ival, flist, glen, header_operators, head, head_l, tail, tail_l, dict)
    pop = sort(pop, by=x -> x.fitness)

    #println("STARTING BEST:")
    #print_indi(pop[1])

    save = pop[1].fitness

    @time while iter < stop && pop[1].fitness > sens

        if iter % 10 == 0 #&& iter > 0
            println(pop[1])
        end

        # Mutation
        m_pop_size = Int(round(pop_size/3))
        mm_pop = Individual[]

        m_pop = gen_m_pop(pop, m_pop_size, 10, "random")
        mt_pop = map(x -> mutate(x, 0.6, 0.2, "change", 0.8), m_pop)
        mm_pop = vcat(mm_pop, mt_pop)

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> mutate(x, 0.6, 0.2, "change", 0.8), m_pop)
        mm_pop = vcat(mm_pop, mt_pop)

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> mutate(x, 0.6, 0.2, "swap", 0.8), m_pop)
        mm_pop = vcat(mm_pop, mt_pop)

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> mutate(x, 0.6, 0.2, "trunc", 0.8), m_pop)
        mm_pop = vcat(mm_pop, mt_pop)

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> mutate(x, 0.6, 0.2, "grow", 0.8), m_pop)
        mm_pop = vcat(mm_pop, mt_pop)

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> mutate(x, 0.6, 0.2, "random", 0.8), m_pop)
        mm_pop = vcat(mm_pop, mt_pop)

        # Mutation of head
        m_pop_size = Int(round(pop_size/3))
        mh_pop = Individual[]

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> muthead(x, 0.6, "jump", 0.8), m_pop)
        mh_pop = vcat(mh_pop, mt_pop)

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> muthead(x, 0.6, "scramble", 0.8), m_pop)
        mh_pop = vcat(mh_pop, mt_pop)

        m_pop = gen_m_pop(pop, m_pop_size, 10, "tournament")
        mt_pop = map(x -> muthead(x, 0.6, "combo", 0.8), m_pop)
        mh_pop = vcat(mh_pop, mt_pop)

        # Crossover
        p_pop_size = Int(round(pop_size/4))
        c_pop = Individual[]

        p_pop = gen_p_pop(pop, p_pop_size, 10, "tournament")
        ct_pop = map(x -> crossover(x..., 0.8, "safe1point", 0.8), p_pop)
        c_pop = vcat(c_pop, ct_pop)

        p_pop = gen_p_pop(pop, p_pop_size, 10, "tournament")
        ct_pop = map(x -> crossover(x..., 0.8, "1point", 0.8), p_pop)
        c_pop = vcat(c_pop, ct_pop)

        p_pop = gen_p_pop(pop, p_pop_size, 10, "tournament")
        ct_pop = map(x -> crossover(x..., 0.8, "2point", 0.8), p_pop)
        c_pop = vcat(c_pop, ct_pop)

        p_pop = gen_p_pop(pop, p_pop_size, 10, "tournament")
        ct_pop = map(x -> crossover(x..., 0.8, "random", 0.8), p_pop)
        c_pop = vcat(c_pop, ct_pop)

        c_pop = vcat(map(x -> x[1], c_pop), map(x -> x[2], c_pop))

        # Creating new pop for next iteration

        #XXX: Doing this really kills diversity very fast.
        pop = vcat(pop, mm_pop, mh_pop, c_pop)

        # Remove entries that are equal according to eq_indi()
        new_pop = Individual[]
        for mem in pop
            if !(true in map(x -> eq_indi(x,mem), new_pop))
                push!(new_pop, mem)
            end
        end
        pop = deepcopy(new_pop)

        # XXX: Make next pop be selected by tournament?

        # XXX: Insert new individuals every now and then?

        # Population for next iteration
        pop = sort(pop, by=x -> x.fitness)
        pop = pop[1:pop_size]

        iter += 1

    end

    println("FINAL BEST (test $i/$len):")
    println(pop[1])

    println("Found in: $iter; time: $(time()-start)")

    push!(iters, iter)
    push!(times, time()-start)

end

println("The times: $times")
println("The iterations: $iters")
println("Time avarage: $(+(times...)/Float64(len))")
println("Iteration avarage: $(+(iters...)/Float64(len))")
