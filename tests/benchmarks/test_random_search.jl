
# Benchmark search: Random search without genetics

include("./../../de_package.jl")

head_l = 5
tail_l = 10
cnum = 4

include("./../../odes/ode1.jl")

pop_size = 400
sens = 10.0^(-7)

times = Float64[]
iters = Int[]
len = 1
for i = 1:len
  iter = 1
  start = time()
  done = false
  best = Inf
  while !done
    pop = gen_pop(1, de, bc, ival, cnum)
    if pop[1].fitness < sens
      done = true
      #print_indi(pop[1])
      #println("FOUND IN: $iter")
    end
    if pop[1].fitness < best
      best = pop[1].fitness
      println(best)
    end
    iter += 1
  end
  push!(iters, iter)
  push!(times, time()-start)
end
println("The times: $times")
println("The iterations: $iters")
println("Time avarage: $(+(times...)/Float64(len))")
println("Iteration avarage: $(+(iters...)/Float64(len))")
