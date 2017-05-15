include("./../../de_package.jl")

head_l = 10
tail_l = 20
cnum = 10

include("./../../odes/ode1.jl")

pop_size = 10

pop = gen_pop(pop_size, de, bc, ival, cnum)

m_pop = map(mutate_head, pop)

for i in 1:length(pop)
  println("-------------------")
  println(pop[i].header)
  println(m_pop[i].header)
  println(pop[i].thestring)
  println(m_pop[i].thestring)
end
