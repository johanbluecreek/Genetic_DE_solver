# Template file for defining a (system of) differential equation(s) with explanations

# Variables
# This is an array of strings containing the names/labels for the variables you want to use.
# Example
vars = ["x", "y"]
# There is a global `VARS` defined that will be used if your own `vars` is not passed to functions explicitly. You may also redefine the global `FLIST`.

# Functions
# This is an array of strings containing the names/labels for the unknown functions you want to solve for.
# Example:
flist = ["y1", "y2"]
# There is a global `FLIST` defined that will be used if your own `flist` is not passed to functions explicitly. You may also redefine the global `FLIST`.

# Differential equations
# This is an array of strings containing the differential equations.
# Syntax for a functions are e.g. "<y1>", with "y1" being the label defined in `flist`.
# Syntax for derivatives are e.g. "<y1xy>", with "y1" being the label in `flist` and `x` and `y` being variables in `vars`.
# Notes on derivatives:
#  * "<y1x>" is the first derviative w.r.t. "x" on "y1", that is "<y1x>" ≡ ∂_x y_1(x, y).
#  * Derivatives are commutative, so even though they will be parsed in a certain order, "<y1xy>" is the same as "<y1yx>"
# Notes on expressions:
#  * Expressions are parsed as strings, meaning if "y1 = x + y" is to be evaluated, "<y1>^2" will result in "x + y^2". So use parentasis, "(<y1>)^2" will result in "(x + y)^2".
#  * The algorithm will look for a solution at which the expression is zero.
# Example:
de = ["(<y1y>)^2 + (<y2>) + sin(x)", "(<y2x>)^2 + (<y1>) + cos(y)"]

# Boundary conditions
# This is an array of arrays of the form `[expression, point]`, where `expression` will be evaluated at `point`. The expression has teh same syntax as for the differential equations, and the points are just an array of floating-numbers following the order of `vars`. The algorithm will look for a solution at which the expression is zero.
# Example:
bc = [["<y1>", [0.0, 1.0]], ["<y2>", [1.0, 0.0]]]

# Interval
# This is an array of tuples of floating-numbers defining the domain of the functions. These are again in the same order as the variables.
# Example:
ival = [(0.0, 0.1), (0.0, 1.0)]

# Global variables
# You may also want to take a look at `src/globaldefault.jl` containing for example the basis set of functions the algorithm will use to try and find a solution.
