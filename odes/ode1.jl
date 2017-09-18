
# Differential equation:
de = ["<ex> - (2 * x - (<e>))/x"]

# Label for the function
flist = ["e"]

# Label for the variables
vars = ["x"]

# Boundary condition
bc = [["(<e>) - 20.1", [("x", 0.1)]]]

# Interval over which to calculate fitness/solve the differential equation
ival = [(0.1,1.0)]

#XXX: These are to be removed:
terminators = vcat(digits, vars, vars, vars, vars, vars);
head = vcat(functions, operators, vars, digits, vars, vars, vars, vars);
tail = vcat(digits, vars);
#XXX: I'm only using them now because of how I wrote testenv.jl and do not want to fix that now...
