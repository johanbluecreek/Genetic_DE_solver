# Benchmark search: Search using all available genetic operators parallelised

###
# Load the package
###
@everywhere include("./../../src/GeneticDESolver.jl")
#XXX: I guess global variables are a really bad thing to use.
# It seems that they have to be loaded by @everywhere to get things to work, which
# is far from convinient.

###
# Load DE to solve
###
@everywhere include("./../../odes/ode1.jl")
#XXX: Gives "UndefVarError: flist not defined"

###
# Set global variables
###

# Basis set
@everywhere functions = ["s", "c", "e", "l", "u"];
@everywhere operators = ["+", "-", "*", "/"];
@everywhere digits = vcat(["$i" for i=range(0,10)], ["p"]);
#XXX: "UndefVarError: functions not defined"
# Relevant variables for PDEs (can also be used as unspecified constants)

# List of the above that terminates an expression tree
@everywhere terminators = vcat(digits, vars, vars, vars, vars, vars);

# header_operators combines genes
@everywhere header_operators = vcat(operators, ["z"]);
# "z" is an operator that deactivates the following chromosome

# Select what should be part of the 'head' of a gene (first head_l entries)
# and which should be part of the 'tail' of a gene (last tail_l entries)
@everywhere head = vcat(functions, operators, vars, digits, vars, vars, vars, vars);
@everywhere tail = vcat(digits, vars);

# Leangth of the head and tail of an elist
@everywhere head_l = 10;
@everywhere tail_l = 20;

# Number of Genes in a Chromosome
@everywhere glen = 4

# Dictionary to translate operators and functions to expressions
@everywhere dict = Dict(
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

        m1_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "change", 0.8), gen_m_pop(pop, m_pop_size, 10, "random"))
        m2_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "change", 0.8), gen_m_pop(pop, m_pop_size, 10, "tournament"))
        m3_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "swap", 0.8), gen_m_pop(pop, m_pop_size, 10, "random"))
        m4_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "trunc", 0.8), gen_m_pop(pop, m_pop_size, 10, "tournament"))
        m5_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "grow", 0.8), gen_m_pop(pop, m_pop_size, 10, "random"))
        m6_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "random", 0.8), gen_m_pop(pop, m_pop_size, 10, "tournament"))

        # Mutation of head
        m_pop_size = Int(round(pop_size/3))

        mh1_pop = @spawn map(x -> muthead(x, 0.6, "jump", 0.8), gen_m_pop(pop, m_pop_size, 10, "tournament"))
        mh2_pop = @spawn map(x -> muthead(x, 0.6, "scramble", 0.8), gen_m_pop(pop, m_pop_size, 10, "tournament"))
        mh3_pop = @spawn map(x -> muthead(x, 0.6, "combo", 0.8), gen_m_pop(pop, m_pop_size, 10, "tournament"))

        # Crossover
        p_pop_size = Int(round(pop_size/4))

        c1_pop = @spawn map(x -> crossover(x..., 0.8, "safe1point", 0.8), gen_p_pop(pop, p_pop_size, 10, "random"))
        c2_pop = @spawn map(x -> crossover(x..., 0.8, "1point", 0.8), gen_p_pop(pop, p_pop_size, 10, "tournament"))
        c3_pop = @spawn map(x -> crossover(x..., 0.8, "2point", 0.8), gen_p_pop(pop, p_pop_size, 10, "tournament"))
        c4_pop = @spawn map(x -> crossover(x..., 0.8, "random", 0.8), gen_p_pop(pop, p_pop_size, 10, "random"))

        # Collecting all new populations
        mm_pop = vcat(fetch(m1_pop), fetch(m2_pop), fetch(m3_pop), fetch(m4_pop), fetch(m5_pop), fetch(m6_pop))

        mh_pop = vcat(fetch(mh1_pop), fetch(mh2_pop), fetch(mh3_pop))

        c_pop = vcat(fetch(c1_pop), fetch(c2_pop), fetch(c3_pop), fetch(c4_pop))
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
