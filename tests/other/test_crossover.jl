
include("./../../de_package.jl")

head_l = 5
tail_l = 10
cnum = 4

include("./../../odes/ode1.jl")

pop_size = 400

pop = gen_pop(pop_size, de, bc, ival, cnum)
c_pop = Individual[]
for i in range(1, 2, Int(round(pop_size)/2))
  c1, c2 = crossover(pop[i], pop[i+1])
  push!(c_pop, c1)
  push!(c_pop, c2)
end

for i in range(1, 2, Int(round(pop_size)/2))
  for j in 1:cnum
    println("parent-pair: $(pop[i].clist[j].clist); $(pop[i+1].clist[j].clist);\n child-pair: $(c_pop[i].clist[j].clist); $(c_pop[i+1].clist[j].clist)")
  end
end
