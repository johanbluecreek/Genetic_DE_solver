# Benchmark search: Search using all available genetic operators parallelised

@everywhere include("./../../de_package.jl")

head_l = 5
tail_l = 10
cnum = 4

include("./../../odes/ode1.jl")

pop_size = 400
stop = 1000
sens = 10.0^(-7)

times = Float64[]
iters = Int[]
len = 30
for i = 1:len

  iter = 1
  start = time()

  pop = gen_pop(pop_size, de, bc, ival, cnum)
  pop = sort(pop, by=x -> x.fitness)

  #println("STARTING BEST:")
  #print_indi(pop[1])

  save = pop[1].fitness

  @time while iter < stop && pop[1].fitness > sens

    if iter % 10 == 0 #&& iter > 0
      println("Top: iteration $iter : ", pop[1].thestring, " : ", pop[1].fitness)
      println("Better: $(save-pop[1].fitness)")
      save = pop[1].fitness
      println("Next:")
      for i in 2:5
        println(pop[i].thestring)
      end
    end

    # Crossover
    p_pop_halfsize = Int(round(pop_size/2))
    #p_pop = gen_p_pop(pop, p_pop_halfsize, "tournament")
    #ct_pop = map(x -> crossover(x..., 0.8, "one-point") , p_pop)
    ct_pop = @spawn map(x -> crossover(x..., 0.8, "one-point") , gen_p_pop(pop, p_pop_halfsize, "tournament"))

    cr_pop = @spawn map(x -> crossover(x..., 0.8, "one-point") , gen_p_pop(pop, p_pop_halfsize, "random"))


    # Mutation
    m_pop_size = Int(round(pop_size/2))
    #m_pop = gen_m_pop(pop, m_pop_size, "random")
    #mc_pop = map(mutate, m_pop)
    mc_pop = @spawn map(mutate, gen_m_pop(pop, m_pop_size, "random"))

    mr_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "random"), gen_m_pop(pop, m_pop_size, "random"))

    mg_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "grow"), gen_m_pop(pop, m_pop_size, "random"))

    mt_pop = @spawn map(x -> mutate(x, 0.6, 0.2, "trunc"), gen_m_pop(pop, m_pop_size, "random"))

    mh_pop = @spawn map(mutate_head, gen_m_pop(pop, m_pop_size, "random"))

    # fetching from Crossover
    nct_pop = Individual[]
    for mem in fetch(ct_pop)
      push!(nct_pop, mem[1])
      push!(nct_pop, mem[2])
    end
    ncr_pop = Individual[]
    for mem in fetch(cr_pop)
      push!(ncr_pop, mem[1])
      push!(ncr_pop, mem[2])
    end

    # combine
    pop = vcat(pop, fetch(mc_pop), fetch(mr_pop), fetch(mg_pop), fetch(mt_pop), fetch(mh_pop), nct_pop, ncr_pop)

    new_pop = Individual[]
    for mem in pop
      if !(true in map(x -> eq_indi(x,mem), new_pop))
        push!(new_pop, mem)
      end
    end
    pop = deepcopy(new_pop)


    # Use tournament to select next pop, to conserve diversity
    #TODO: Drop the selected one from pop instead, so that this does not loop so much
    new_pop = Individual[]
    pool = deepcopy(pop)
    while length(new_pop) < pop_size
      indi = m_select(pool, "tournament")
      filter!( e -> !(eq_indi(e, indi)), pool)
      push!(new_pop, indi)
    end
    pop = deepcopy(new_pop)

    pop = sort(pop, by=x -> x.fitness)

    # pop = pop[1:pop_size]

    # To keep diversity up, insert new Individuals
    if iter % 3 == 0
      cut = Int(round(pop_size/4))
      pop = vcat(pop[1:(pop_size - cut)], gen_pop(cut, de, bc, ival, cnum))
    end

    iter += 1

  end

  #println("FINAL BEST:")
  #print_indi(pop[1])

  println("Found in: $iter; time: $(time()-start)")

  push!(iters, iter)
  push!(times, time()-start)

end

println("The times: $times")
println("The iterations: $iters")
println("Time avarage: $(+(times...)/Float64(len))")
println("Iteration avarage: $(+(iters...)/Float64(len))")
