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
m_pop = map(x -> mutate(x, 1.0, 1.0, "swap"), pop)


for i in 1:pop_size
  println(pop[i].thestring)
  println(m_pop[i].thestring)
end
