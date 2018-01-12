
de = "((<e>)*((<exx>)+(<eyy>)) - (<ex>)^2 - (<ey>)^2)/(<e>)^2 - (<e>)^2"

# Label for the function
flist = ["e"]

# Label for the variables
vars = ["x", "y"]

bc = [["0", [("x", 0.0), ("y", 0.0)]]]

ival = [(0.1, 100.0), (0.1, 100.0)]
