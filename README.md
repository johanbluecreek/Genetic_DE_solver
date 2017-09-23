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

## Overview

The solutions are found using the following structure

<table>
  <tr>
    <td colspan="13"><b>Individual</b><br/> Represents solution to system of differential equations </td>
  </tr>
  <tr>
    <td colspan="4"><b>Chromosome</b><br/> Expression for Function <code>1</code> </td>
    <td colspan="4"><b>Chromosome</b><br/> Expression for Function <code>2</code> </td>
    <td ><b>...</b></td>
    <td colspan="4"><b>Chromosome</b><br/> Expression for Function <code>n</code></td>
  </tr>
  <tr>
    <td><b>Gene</b><br/> Sub-expression <code>1</code> </td>
    <td><b>Gene</b><br/> Sub-expression <code>2</code> </td>
    <td><b>...</b></td>
    <td><b>Gene</b><br/> Sub-expression <code>m</code> </td>
    <td colspan="4"><b>Genes...</b></td>
    <td><b>...</b></td>
    <td colspan="4"><b>Genes...</b></td>
  </tr>
</table>

That is, the smallest components are the mathematical expressions that the `Gene`-types give rise to. These genes are combined in the `Chromosome`-types to form a complete mathematical expression representing a function. The `Individual`-type consists of all the chromosomes and represents a complete solution to a system of differential equations.

The above table depicts the structure of a solution to a system of `n` differential equations/functions, where each function is made up by `m` sub-expressions.

### Mutation

There are several mutation operators implemented, and they are divided in two parts.

**Gene mutation:** The main form of mutations acts on the `Gene`s and changes the mathematical expressions they give rise to. See `?mutate` for further explanation on this.

**Chromosome mutation:** A supplementary form of mutation that changes how the `Gene`s are combined to form a `Chromosome`, leaving the mathematical expressions of the `Gene`s fixed. See `?muthead` for more information.

### Crossover

There are some crossover operators implemented. These are interactions between `Individual`s such that they produce new solutions as offspring. See `?crossover` for more information.

### Fitness calculation

The fitness if calculated similar to that of **[TL]**. The fitness is here divided in three parts

* **error:** Sum of squares of the differential equations evaluated at a set number of uniformly distributed points over a given interval. (Also in **[TL]**)
* **penalty:** Absolute value of the mismatch for the boundary condition. (Also in **[TL]**)
* **shape:** This is the same as error, but evaluated for derivatives of the differential equation, with a decay-factor for each derivative. (Not in **[TL]**)

both penalty and shape can be turned off at your preference. The fitness is then simply given as the sum of these parts.

### Comments

The system-number `n` is actually only set by the number of functions involved in the system, meaning that the solver is blind to under- or over-determined systems.

That is, in an under-determined system of, say, five functions and one differential equation, the result will be *a* solution (if one exists, else an approximation) to this "system". In an over-determined system of, say, one function and five differential equations, the result will be the best approximation (or possibly exact solution iff the system is degenerate) found when the termination conditions were met.

The sub-expression number `m` is independent of any details of the problem, and can be chosen at will. It should however be noted that larger `m` will have negative effect on the time it takes to execute each iteration. See **[MF]** regarding the introduction of sub-expressions.

## Usage

At the present time, the documentation is the code, unfortunately. But having Julia installed properly, one can execute the files under `tests/` using `julia` to try it out. This will also tell you what packages you are missing (that is, requirements), but this should be it:

 * Calculus
 * Iterators

install with running `Pkg.add()` in `julia`.

## Contributing

Any contribution would be appreciated, you can start off small, e.g.

* Write the files for any of the differential equation in **[TL]**.

try it out and tune it, e.g.

* Find out optimal settings for mutation (rate, selection, sizes, operations), crossover, etc. for the different differential equations, for example `ode3.jl`.

or help developing, for example

* See the TODO below
* Read TODO, XXX, and FIXME comments in the code

## Repository content

* `odes/`: Folder for ODE:s from **[TL]**.
* `pdes/`: Folder for PDE:s from **[TL]** and others.
* `tests/`: Folder for testing various aspects of the code.
* `src/`: Source files.
* `src/GeneticDESolver.jl`: Main file that makes all functions available. Type `include("src/GeneticDESolver.jl")` to load.
* `src/genops/`: Genetic operators.
* `src/testenv.jl`: Small test environment to test out different features. Load with `include("src/GeneticDESolver.jl"); include("src/testenv.jl")`.

## TODO/Known issues

Crude TODO-list/feature wish-list

* Implement the rest of the crossover operators listed in **[MF]** and **[E]**
* Also implement their parent-selection methods
* Add the rest of the differential equations in **[TL]**
* Document functions and code
* Find optimal settings for solving the differential equations of **[TL]**.
* Make a pedagogical introduction in a jupyter-notebook for how all this works.

Known issues

* There is a problem, visible in for example test_everything_parallel.jl, that it takes longer and longer time to complete one iteration as it is progressing. There is also memory leaks for this file.
