
include("./../../src/GeneticDESolver.jl")
include("./../../odes/ode1.jl")

# Rust voodoo
ln(x) = log(x)

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
    "l" => "ln(<expr>)",
    "e" => "exp(<expr>)",
    "u" => "(<expr>)",
);
const pop_size = 10

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
stop = 2000
tworst = 0.0
times = Float64[]
bworst = 0
bytess = []
gctimes = []

tavgs = []

while iter < stop

    val, tsave, bytes, gctime, memallocs = @timed npop = map(x -> mutate(x, 0.05, 1.0, "random", 1.0), pop)

    push!(times, tsave)
    push!(bytess, bytes)
    push!(gctimes, gctime)


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
