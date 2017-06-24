# Genetic DE solver

This is a differential equations solver that finds analytical solutions using genetic programming. This is written in [Julia](https://github.com/JuliaLang/julia) and relies on the [Calculus](https://github.com/johnmyleswhite/Calculus.jl) package for symbolic derivatives.

Compared to [GP_DE_solver (development halted)](https://github.com/johanbluecreek/GP_DE_solver) written in Python, this is faster, more reliable, have more genetic operators implemented, and uses a parse-tree sort of structure rather than a grammar.

## Introduction

The inspiration for this project is taken from

* **[TL]** TSOULOS, Ioannis G. et LAGARIS, Isaac E. "*Solving differential equations with genetic programming*". Genetic Programming and Evolvable Machines, 2006, vol. 7, no 1, p. 33-54.

the following books have been used as a reference

* **[MF]** "*How to Solve It: Modern Heuristics*" by Zbigniew Michalewicz, David
B. Fogel
* **[E]** "*Computational Intelligence: An Introduction*" by Andries P.
Engelbrecht

## Usage

At the present time, the documentation is the code, unfortunately. But having Julia installed properly, one can execute the files under `tests/` using `julia` to try it out.

## Contributing

Any contribution would be appreciated, you can start off small, e.g.

* Write the files for any of the differential equation in **[TL]** (note, systems and PDEs have not yet been fully implemented, hence syntax not yet fixed)

try it out and tune it, e.g.

* Find out optimal settings for mutation (rate, selection, sizes, operations), crossover, etc. for the different differential equations, for example ode3.jl.

or help developing, for example

* See the TODO below
* Read TODO, XXX, and FIXME comments in the code
* Highest priority right now is implmenting support for systems of differential equations. Follow development in branch 'system'

## Repository content

* `odes/`: Folder for ODE:s from **[TL]**.
* `tests/`: Folder for testing various aspects of the code.
* `de_package.jl`: File containing the types, functions, and genetic operators.

## TODO/Known issues

Crude TODO-list/feature wish-list

* Implement the rest of the crossover operators listed in **[MF]** and **[E]**
* Also implement their parent-selection methods
* Add support for systems of differential equations (see branch 'system')
* Add support for PDE:s (see branch 'partial', and 'system')
* Add the rest of the differential equations in **[TL]**
* Document functions and code

Known issues

* It still does not solve all ODE:s (e.g. `ode3`) in reasonable time (if at all).
