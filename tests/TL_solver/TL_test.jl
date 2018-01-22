
include("../../solvers/TL_solver.jl")

diffs = ["../../des/odes/ode1.jl"]
sols = []

loops = 5

for diff in diffs
    for i in 1:loops
        sol = TL_solver(diff, 100, Inf, 1000, 100, 10.0^(-7))
        append!(sols, [sol])
    end
end

exprs = map(x->x[1], sols)
times = map(x->x[2], sols)
iters = map(x->x[3], sols)

fails = map(x->isequal(x, Inf), exprs)

deleteat!(times, fails)
deleteat!(iters, fails)

println(" ")
println("Sucesses: ", length(times))
println("Min time: ", min(times...))
println("Max time: ", max(times...))
println("Avg time: ", (sum(times)/length(times)))
println("Min iter: ", min(iters...))
println("Max iter: ", max(iters...))
println("Avg iter: ", (sum(iters)/length(iters)))
