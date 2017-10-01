
# XXX: After some research:
#=
This code seems to be affected by:
https://github.com/JuliaLang/julia/issues/19013
https://github.com/JuliaLang/julia/issues/18446
https://github.com/JuliaLang/julia/issues/14495
... and possily refs within.

=#

include("./../../src/GeneticDESolver.jl")
include("./../../odes/ode1.jl")

const functions = ["s", "c", "e", "l"];
const operators = ["+", "-", "*", "/"];
const digits = vcat(["$i" for i=range(0,10)]);
const terminators = vcat(digits, vars);
const header_operators = vcat(operators);
const head = vcat(repeat(operators, outer=5), repeat(["u"], outer=20), repeat(functions, outer=5), repeat(digits, outer=2), repeat(vars, outer=20));
const tail = vcat(digits, repeat(vars, outer=10));
const head_l = 10;
const tail_l = head_l*2;
const glen = 2
const dict = Dict(
    "+" => "(<expr>)+(<expr>)",
    "-" => "(<expr>)-(<expr>)",
    "*" => "(<expr>)*(<expr>)",
    "/" => "(<expr>)/(<expr>)",
    "s" => "sin(<expr>)",
    "c" => "cos(<expr>)",
    "l" => "log(<expr>)",
    "e" => "exp(<expr>)",
    "u" => "(<expr>)",
);
const pop_size = 10

#XXX: The below does not show the same behaviour
#XXX: But why did it, when I exchanged fisk/defunc ?
#XXX: Because it somehow changes the func in all Individuals?
#XXX: Will simply removing the defunc functionality solve it?
#XXX: THIS IS FOR TOMORROW! THATS TODAY TO YOU!
function mokgenind()
    # gen_indi generates a function, so lets do that.
    y = rand(1:10)
    f_body = :(x^$y)
    f_call = Expr(:call,:fisk,:x)
    eval(Expr(:function,f_call,f_body))
    a = Base.invokelatest(fisk, 2)

    return fisk, a
end

function mokrepind(xin)
    # reparse_indi makes a copy and calls gen_indi, so do that.
    x = deepcopy(xin)
    f, a = mokgenind()
    return x
end

function mokmut(xin)
    # mutation makes a copy and calls reparse_indi
    x = deepcopy(xin)
    x = mokrepind(x)
    return x
end

pop = gen_pop(pop_size, de, bc, ival, flist, glen, header_operators, head, head_l, tail, tail_l, dict)

iter = 1
stop = 2001
tworst = 0.0
times = Float64[]
bworst = 0
bytess = []
gctimes = []

tavgs = []

while iter < stop

    val, tsave, bytes, gctime, memallocs = @timed npop = map(x -> mutate(x, 0.05, 1.0, "random", 1.0), pop)

    # val, tsave, bytes, gctime, memallocs = @timed npop = mokmut(pop)

    push!(times, tsave)
    push!(bytess, bytes)
    push!(gctimes, gctime)

#=
    if tsave > tworst && iter > 1
        tworst = tsave
        println("$iter had worse time: $tworst")
    end

    if bytes > bworst && iter > 1
        bworst = bytes
        println("$iter used more memory: $bworst")
    end

    if iter % 100 == 0
        println("avg time: $(sum(times)/length(times))")
        println("avg bytes: $(sum(bytess)/length(bytess))")
        println("avg gctime: $(sum(gctimes)/length(gctimes))")
        times = Float64[]
        bytess = []
    end

=#

if iter % 100 == 0
    tavg = (sum(times)/length(times))
    if !(false in (tavg .> tavgs))
        println("$iter: worse!")
    end
    push!(tavgs, tavg)
    times = Float64[]
end

    iter += 1
end
#XXX: Looping npop = map(x -> mutate(x, 0.05, 1.0, "random", 1.0), pop)
# gets worse by a factor ~2 per 2000 iterations (in time).
# 'bytes' seems to be stable
# Looking at htop during execution shows a steady increase in memory
#=
    $ julia memoryloops.jl
2 was worse: 0.32207298278808594
3 was worse: 0.39785194396972656
avg: 0.29034024715423584
avg: 0.2721996736526489
avg: 0.2882564735412598
369 was worse: 0.4127919673919678
avg: 0.3124001216888428
455 was worse: 0.4817500114440918
avg: 0.3345027780532837
avg: 0.3529787945747376
avg: 0.38689337491989134
728 was worse: 0.4966590404510498
avg: 0.417019464969635
avg: 0.42663576602935793
924 was worse: 0.5299191474914551
943 was worse: 0.5307579040527344
avg: 0.4614319944381714
1031 was worse: 0.532965898513794
1033 was worse: 0.555980920791626
1056 was worse: 0.6265590190887451
avg: 0.4788250732421875
avg: 0.5151924085617066
1218 was worse: 0.7778811454772949
avg: 0.5427266073226928
avg: 0.555482611656189
avg: 0.5702780961990357
avg: 0.5973370409011841
avg: 0.6241191744804382
avg: 0.6455622100830078
1841 was worse: 0.9063899517059326
avg: 0.6748245382308959
1956 was worse: 0.906782865524292
=#
#=
    $ julia memoryloops.jl
2 had worse time: 0.352611984
2 used more memory: 11314625
3 had worse time: 0.37092432
3 used more memory: 11546281
4 had worse time: 0.381921997
4 used more memory: 11951764
27 had worse time: 0.40270329
27 used more memory: 12364540
avg time: 0.30566370014
avg bytes: 1.048699798e7
avg time: 0.29869712839
avg bytes: 9.90928983e6
avg time: 0.33627694853000006
avg bytes: 9.89405258e6
316 had worse time: 0.412273114
350 had worse time: 0.414259717
359 had worse time: 0.429652135
avg time: 0.368413629
avg bytes: 9.98067597e6
414 had worse time: 0.508133567
497 had worse time: 0.577097192
avg time: 0.37241267445000004
avg bytes: 9.96209988e6
avg time: 0.40187657778
avg bytes: 9.94562568e6
avg time: 0.4356459253600001
avg bytes: 9.96627483e6
738 had worse time: 0.588692506
753 had worse time: 0.663094578
avg time: 0.46983433046
avg bytes: 9.88800332e6
avg time: 0.48688681631999997
avg bytes: 9.92168447e6
avg time: 0.50081427905
avg bytes: 1.003341807e7
1084 had worse time: 0.715969249
avg time: 0.5257692937
avg bytes: 9.8870916e6
avg time: 0.5502780686400001
avg bytes: 9.97378732e6
1201 had worse time: 0.759891901
1226 had worse time: 0.771579676
avg time: 0.58003076007
avg bytes: 9.98452123e6
1356 had worse time: 0.785833847
1380 had worse time: 0.794682189
avg time: 0.61212571816
avg bytes: 9.94223305e6
avg time: 0.61047666049
avg bytes: 9.87148341e6
1520 had worse time: 0.87408254
avg time: 0.6582714728600001
avg bytes: 9.92515136e6
avg time: 0.6875104116700002
avg bytes: 9.94537203e6
1716 had worse time: 0.875687277
1771 had worse time: 0.938985053
avg time: 0.7328375529900001
avg bytes: 1.000354221e7
1891 had worse time: 0.948846164
avg time: 0.75496086376
avg bytes: 9.94637072e6
1904 had worse time: 0.964578244
1906 had worse time: 0.969331025
1907 had worse time: 1.0266973
1914 had worse time: 1.07298249
1918 had worse time: 1.07493805
1920 had worse time: 1.277098666
1924 had worse time: 1.374999605

note that at ~1900 I started compiling code on my computer and ruined the test. Still, increasing time up to that point.
=#

#XXX Nothing seem to have changed by declaring stuff as `const`

#XXX: IT also keeps hanging... why. Oh, old calculus? yes. Pkg.checkout("Calculus") solves it.

#XXX: gctime also increasing, but cannot account for everything.

#XXX: can I make a pretend function that mimics this? Say I...

#XXX: it is already here:
#=
$ julia memoryloops.jl
100: worse!
800: worse!
900: worse!
1000: worse!
1100: worse!
1200: worse!
1300: worse!
1400: worse!
1500: worse!
1600: worse!
1700: worse!
1800: worse!
1900: worse!
2000: worse!
=#
#XXX: Only doing deep copies of a pop.

#XXX: I removed the use of defunc, so that is nolonger around, but memory leak is still there and slowdown is... still there.

#XXX: Make defunc anonamyous
#XXX: This does not help
