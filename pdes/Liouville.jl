# Gives solutions at contant -Inf
#de = "(<exx>) + (<eyy>) - exp(2*(<e>))"

# Zero is a solution here...
# <e> -> log(<e>)
de = "((<e>)*((<exx>)+(<eyy>)) - (<ex>)^2 - (<ey>)^2)/(<e>)^2 - (<e>)^2"

bc = ["0", 0.0]
ival = [0.0, 100.0]
