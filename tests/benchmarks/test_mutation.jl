include("./../../de_package.jl")

head_l = 5
tail_l = 10
cnum = 4

include("./../../odes/ode1.jl")

pop_size = 20

pop = gen_pop(pop_size, de, bc, ival, cnum)

# Simply change the method-string below and run program.
# The original versus mutated population will be displayed pairwise for you
# to compare.

#m_pop = map(x -> mutate(x, 1.0, 1.0, "swap"), pop)
m_pop = map(mutate_jump, pop)

for i in 1:pop_size
  println("----")
  println("head: ", pop[i].header)
  println("ori: ", pop[i].thestring)
  println("new: ", m_pop[i].thestring)
end
