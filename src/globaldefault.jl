# Global Variables
# These are the variables that will be used by functions that rely on them if nothing else is passed to those functions.

#TODO: Be more strict with the types

FUNCTIONS = ["s", "c", "e", "l", "u"];
OPERATORS = ["+", "-", "*", "/"];
DIGITS = vcat(["$i" for i=range(0,10)]);

HEAD = [FUNCTIONS; OPERATORS; DIGITS];
TAIL = DIGITS;

HEAD_L = 5;
TAIL_L = 6;

HEADER_OPERATORS = [OPERATORS; ["z"]];

PENALTY_FACTOR = 100;
SHAPE_DECAY = 1000;

DICT = Dict(
    "+" => "(<expr>)+(<expr>)",
    "-" => "(<expr>)-(<expr>)",
    "*" => "(<expr>)*(<expr>)",
    "/" => "(<expr>)/(<expr>)",
    "s" => "sin(<expr>)",
    "c" => "cos(<expr>)",
    "l" => "log(<expr>)",
    "e" => "exp(<expr>)",
    "u" => "(<expr>)",
);

FLIST = ["e"];
VARS = ["x", "y"];

TERMINATORS = [DIGITS; VARS];

GLEN = 2;
