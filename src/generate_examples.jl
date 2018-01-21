# Generate examples for function documentation

include("./GeneticDESolver.jl")

function print_exec(exec)
    println("julia> $(exec)")
    println(eval(parse(exec)))
    println("")
end

println("")
println("<<<<<<<<<<<<<<<<<< genfunc.jl >>>>>>>>>>>>>>>>>>")
println("")

println("----------- safe_string -----------")

exec = "safe_string(\"2x\")"
print_exec(exec)

println("----------- init_elist -----------")

exec = "init_elist()"
print_exec(exec)

println("----------- parse_elist -----------")

exec = "elist = init_elist()"
print_exec(exec)

exec = "parse_elist(elist)"
print_exec(exec)

println("----------- parse_tree -----------")

exec = "elist = init_elist()"
print_exec(exec)

exec = "parse_tree(elist)"
print_exec(exec)

println("----------- init_gene -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "typeof(g)"
print_exec(exec)

println("----------- reparse_gene -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "g.elist = init_elist()"
print_exec(exec)

exec = "reparse_gene(g)"
print_exec(exec)

println("----------- init_chromo -----------")

exec = "c = init_chromo()"
print_exec(exec)

exec = "c.glist"
print_exec(exec)

println("----------- reparse_chromo -----------")

exec = "c = init_chromo()"
print_exec(exec)

exec = "c.glist = [init_gene(), init_gene()]"
print_exec(exec)

exec = "reparse_chromo(c)"
print_exec(exec)

println("----------- eq_chromo -----------")

exec = "c1 = init_chromo()"
print_exec(exec)

exec = "c2 = init_chromo()"
print_exec(exec)

exec = "eq_chromo(c1,c2)"
print_exec(exec)

exec = "eq_chromo(c1,c1)"
print_exec(exec)

println("----------- gen_indi_attributes -----------")

exec = "c1 = init_chromo()"
print_exec(exec)

exec = "c2 = init_chromo()"
print_exec(exec)

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "gen_indi_attributes([c1, c2], de, bc, ival)"
print_exec(exec)

println("----------- init_indi -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "i = init_indi(de, bc, ival)"
print_exec(exec)

exec = "typeof(i)"
print_exec(exec)

println("----------- reparse_indi -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "i = init_indi(de, bc, ival)"
print_exec(exec)

exec = "i.clist[1] = init_chromo()"
print_exec(exec)

exec = "reparse_indi(i)"
print_exec(exec)

println("----------- parse_expr -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "i = init_indi(de, bc, ival)"
print_exec(exec)

exec = "parse_expr(de, i)"
print_exec(exec)

println("")
println("<<<<<<<<<<<<<<<<<< types.jl >>>>>>>>>>>>>>>>>>")
println("")

println("----------- Gene -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "typeof(g)"
print_exec(exec)

exec = "g.elist"
print_exec(exec)

exec = "g.thestring"
print_exec(exec)

exec = "g.tree"
print_exec(exec)

println("----------- Chromosome -----------")

exec = "c = init_chromo()"
print_exec(exec)

exec = "typeof(c)"
print_exec(exec)

exec = "c.glist"
print_exec(exec)

exec = "c.header"
print_exec(exec)

exec = "c.thestring"
print_exec(exec)

println("----------- Individual -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "i = init_indi(de, bc, ival)"
print_exec(exec)

exec = "typeof(i)"
print_exec(exec)

exec = "i.clist"
print_exec(exec)

exec = "i.error"
print_exec(exec)

exec = "i.penalty"
print_exec(exec)

exec = "i.shape"
print_exec(exec)

exec = "i.fitness"
print_exec(exec)

println("")
println("<<<<<<<<<<<<<<<<<< genops/crossover.jl >>>>>>>>>>>>>>>>>>")
println("")

println("----------- cross_safe1point -----------")

exec = "c1 = init_chromo()"
print_exec(exec)

exec = "c2 = init_chromo()"
print_exec(exec)

exec = "c3, c4 = cross_safe1point(c1, c2, 1.0)"
print_exec(exec)

exec = "[map(x -> x.elist, c1.glist), map(x -> x.elist, c2.glist)]"
print_exec(exec)

exec = "[map(x -> x.elist, c3.glist), map(x -> x.elist, c4.glist)]"
print_exec(exec)

println("----------- cross_1point -----------")

exec = "c1 = init_chromo()"
print_exec(exec)

exec = "c2 = init_chromo()"
print_exec(exec)

exec = "c3, c4 = cross_1point(c1, c2, 1.0)"
print_exec(exec)

exec = "[map(x -> x.elist, c1.glist), map(x -> x.elist, c2.glist)]"
print_exec(exec)

exec = "[map(x -> x.elist, c3.glist), map(x -> x.elist, c4.glist)]"
print_exec(exec)

println("----------- cross_2point -----------")

exec = "c1 = init_chromo()"
print_exec(exec)

exec = "c2 = init_chromo()"
print_exec(exec)

exec = "c3, c4 = cross_2point(c1, c2, 1.0)"
print_exec(exec)

exec = "[map(x -> x.elist, c1.glist), map(x -> x.elist, c2.glist)]"
print_exec(exec)

exec = "[map(x -> x.elist, c3.glist), map(x -> x.elist, c4.glist)]"
print_exec(exec)

println("----------- cross_random -----------")

exec = "glen = 5"
print_exec(exec)

exec = "c1 = init_chromo(glen)"
print_exec(exec)

exec = "c2 = init_chromo(glen)"
print_exec(exec)

exec = "c3, c4 = cross_random(c1, c2, 0.5)"
print_exec(exec)

exec = "[map(x -> x.elist, c1.glist), map(x -> x.elist, c2.glist)]"
print_exec(exec)

exec = "[map(x -> x.elist, c3.glist), map(x -> x.elist, c4.glist)]"
print_exec(exec)

println("----------- crossover -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "i1 = init_indi(de, bc, ival)"
print_exec(exec)

exec = "i2 = init_indi(de, bc, ival)"
print_exec(exec)

exec = "i3, i4 = crossover(i1, i2, 1.0, \"safe1point\", 1.0)"
print_exec(exec)

println("")
println("<<<<<<<<<<<<<<<<<< genops/mutation.jl >>>>>>>>>>>>>>>>>>")
println("")

println("----------- mut_change -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "mut_change(g, 1.0)"
print_exec(exec)

println("----------- mut_random -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "mut_random(g, 1.0)"
print_exec(exec)

println("----------- mut_trunc -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "mut_trunc(g)"
print_exec(exec)

println("----------- mut_grow -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "mut_grow(g)"
print_exec(exec)

println("----------- mut_swap -----------")

exec = "g = init_gene()"
print_exec(exec)

exec = "mut_swap(g)"
print_exec(exec)

println("----------- mutation -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "i = init_indi(de, bc, ival)"
print_exec(exec)

exec = "mutate(i, 0.6, 0.2, \"change\", 0.8)"
print_exec(exec)

println("")
println("<<<<<<<<<<<<<<<<<< genops/muthead.jl >>>>>>>>>>>>>>>>>>")
println("")

println("----------- muth_scramble -----------")

exec = "glen = 5; c = init_chromo(glen)"
print_exec(exec)

exec = "muth_scramble(c)"
print_exec(exec)

println("----------- muth_jump -----------")

exec = "glen = 5; c = init_chromo(glen)"
print_exec(exec)

exec = "muth_jump(c)"
print_exec(exec)

println("----------- muth_combo -----------")

exec = "glen = 5; c = init_chromo(glen)"
print_exec(exec)

exec = "muth_combo(c)"
print_exec(exec)

println("----------- muthead -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "i = init_indi(de, bc, ival)"
print_exec(exec)

exec = "muthead(i, 1.0, \"combo\", 1.0)"
print_exec(exec)

println("")
println("<<<<<<<<<<<<<<<<<< genops/population.jl >>>>>>>>>>>>>>>>>>")
println("")

println("----------- gen_pop -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "gen_pop(5, de, bc, ival)"
print_exec(exec)

println("")
println("<<<<<<<<<<<<<<<<<< genops/parentsel.jl >>>>>>>>>>>>>>>>>>")
println("")

println("----------- pmet_tournament -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "pop = gen_pop(100, de, bc, ival)"
print_exec(exec)

exec = "pmet_tournament(pop, 5)"
print_exec(exec)

println("----------- pmet_random -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "pop = gen_pop(100, de, bc, ival)"
print_exec(exec)

exec = "pmet_random(pop)"
print_exec(exec)

println("----------- p_select -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "pop = gen_pop(100, de, bc, ival)"
print_exec(exec)

exec = "p_select(pop, 5, \"tournament\")"
print_exec(exec)

println("----------- gen_p_pop -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "pop = gen_pop(100, de, bc, ival)"
print_exec(exec)

exec = "gen_p_pop(pop, 10, 5, \"tournament\")"
print_exec(exec)

println("----------- gen_m_pop -----------")

exec = "include(\"../des/de_template.jl\"); FLIST = flist; [de, bc, ival]"
print_exec(exec)

exec = "pop = gen_pop(100, de, bc, ival)"
print_exec(exec)

exec = "gen_m_pop(pop, 10, 5, \"tournament\")"
print_exec(exec)
