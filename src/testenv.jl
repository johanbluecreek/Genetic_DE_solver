
#=
This is here to provide a test enviornment
=#

# Basis set
functions = ["s", "c", "e", "l", "u"]
operators = ["+", "-", "*", "/"]
digits = vcat(["$i" for i=range(0,10)], ["p"])
# Relevant variables for PDEs (can also be used as unspecified constants)
vars = ["x", "y"]
# Relevant function labels for systems
flist = ["f1", "f2"]

terminators = vcat(digits, vars, vars, vars, vars, vars)

# "z" is an operator that deactivates the following chromosome
header_operators = vcat(operators, ["z"])

head = vcat(functions, operators, vars, digits, vars, vars, vars, vars)
tail = vcat(digits, vars)

head_l = 10
tail_l = 20

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
)
