
de = ["(<exx>) + (<eyy>) - exp(-x)*(x-2+y^3+6x)"]

# Label for the function
flist = ["e"]

# Label for the variables
vars = ["x", "y"]

bc = [
    ["(<e>)", [0.0, "y"]],
    ["(<e>) - sin(1)*cos(y)", [1.0, "y"]],
    ["(<e>) - sin(x)", ["x", 0.0]],
    ["(<e>) - sin(x)*cos(1)", ["x", 1.0]]
]

ival = [(0.0, 0.1), (0.0, 0.1)]
