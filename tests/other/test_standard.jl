include("./../../de_package.jl")

head_l = 10
tail_l = 20
cnum = 10

include("./../../pdes/pde1.jl")

pop_size = 10

stop = 20
sens = 10.0^(-10)

pop = gen_pop(pop_size, de, bc, ival, cnum)

for mem in pop
  println("--------------------")
  println(mem.thestring)
  println(mem.fitness)
end
println("-------------------")


