# Differential equation:
de = ["cos(x) + (<y1>)^2 + (<y2>) - (x^2 + sin(x)^2) - (<y1x>)", "2*x - x^2 * sin(x) + (<y1>)*(<y2>) - (<y2x>)"]

# Label for the function
flist = ["y1", "y2"]

# Label for the variables
vars = ["x"]

# Boundary condition
bc = [["<y1>", [0.0]], ["<y2>", [0.0]]]

# Interval over which to calculate fitness/solve the differential equation
ival = [(0.0,1.0)]
