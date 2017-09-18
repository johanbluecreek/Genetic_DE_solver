
#=
This is here to provide a test enviornment
=#

# Basis set
functions = ["s", "c", "e", "l", "u"];
operators = ["+", "-", "*", "/"];
digits = vcat(["$i" for i=range(0,10)], ["p"]);
# Relevant variables for PDEs (can also be used as unspecified constants)
vars = ["x", "y"];
# Relevant function labels for systems
flist = ["f1", "f2"];

# List of the above that terminates an expression tree
terminators = vcat(digits, vars, vars, vars, vars, vars);

# header_operators combines genes
header_operators = vcat(operators, ["z"]);
# "z" is an operator that deactivates the following chromosome

# Select what should be part of the 'head' of a gene (first head_l entries)
# and which should be part of the 'tail' of a gene (last tail_l entries)
head = vcat(functions, operators, vars, digits, vars, vars, vars, vars);
tail = vcat(digits, vars);

head_l = 10;
tail_l = 20;

# Dictionary to translate operators and functions to expressions
dict = Dict(
    "+" => "(<expr>)+(<expr>)",
    "-" => "(<expr>)-(<expr>)",
    "*" => "(<expr>)*(<expr>)",
    "/" => "(<expr>)/(<expr>)",
    "s" => "sin(<expr>)",
    "c" => "cos(<expr>)",
    "l" => "log(<expr>)",
    "e" => "exp(<expr>)",
    "u" => "(<expr>)",
    "p" => "pi"
);

# Below we generate the necessary variables to run tests and examples in the docs
# These require all functions, that is, the package to be run. If it has, but these
# are undesireable, set GDES_loaded=false
if isdefined(:GDES_loaded)
    if GDES_loaded
        # These are do not require the package functions, but are also only useful in a
        # strict testing env.
        test_expr = "-(<e>)+(<exx>)^2+(<eyy>)^2"
        test_exprs = ["-(<f1>)+(<f2>)+(<f1xx>)^2+(<f1yy>)^2", "(<f1>)-(<f2>)+(<f2xx>)^2+(<f2yy>)^2"]
        test_de = test_exprs
        test_bc = [
            ["(<f1x>) - y", [("x", 0.1), ("y", 0.2)]],
            ["(<f1y>) - x", [("x", 0.2), ("y", 0.3)]]
        ]
        test_ival = [(0.1,1.0), (0.2,1.1)]

        test_elist = init_elist()
        test_gene = init_gene()
        test_chromo = init_chromo()
        test_chromo1, test_chromo2 = init_chromo(), init_chromo()
        test_indi = init_indi(test_de, test_bc, test_ival)
        test_indi1, test_indi2 = init_indi(test_de, test_bc, test_ival), init_indi(test_de, test_bc, test_ival)

        test_pop = gen_pop(5, test_de, test_bc, test_ival)
    else
        println("Info: Only defining most basic test enviornment.")
    end
else
    println("Warn: GeneticDeSolver.jl appears to not have been included. Only basic test enviornment will be present.")
end
