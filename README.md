# Genetic DE solver

This is a differential equations solver (ODE only, presently) that finds analytical solutions using genetic programming. This is written in [Julia](https://github.com/JuliaLang/julia) and relies on the [Calculus](https://github.com/johnmyleswhite/Calculus.jl) package for symbolic derivatives.

Compared to [GP_DE_solver](https://github.com/johanbluecreek/GP_DE_solver) written in Python, this is faster, more reliable, have more genetic operators implemented, and uses a parse-tree sort of structure rather than a grammar.

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
    <td colspan="13">**Individual**<br/> Represents solution to system of differential equations </td>
  </tr>
  <tr>
    <td colspan="4">**Chromosome**<br/> Expression for Function `1` </td>
    <td colspan="4">**Chromosome**<br/> Expression for Function `2` </td>
    <td >**...**</td>
    <td colspan="4">**Chromosome**<br/> Expression for Function `n`</td>
  </tr>
  <tr>
    <td>**Gene**<br/> Sub-expression `1` </td>
    <td>**Gene**<br/> Sub-expression `2` </td>
    <td>**...**</td>
    <td>**Gene**<br/> Sub-expression `m` </td>
    <td colspan="4">**Genes...**</td>
    <td>**...**</td>
    <td colspan="4">**Genes...**</td>
  </tr>
</table>

That is, the smallest components are the mathematical expressions that the `Gene`-types give rise to. These genes are combined in the `Chromosome`-types to form a complete mathematical expression representing a function. The `Individual`-type consists of all the chromosomes and represents a complete solution to a system of differential equations.

The above table depicts the structure of a solution to a system of `n` differential equations/functions, where each function is made up by `m` sub-expressions.

### Mutation

There are several mutation operators implemented, and they are divided in two parts.

**Gene mutation:** The main form of mutations acts on the `Gene`s and changes the mathematical expressions they give rise to. See `?mutate` for further explanation on this.

**Chromosome mutation:** A supplementary form of mutation that changes how the `Gene`s are combined to form a `Chromosome`, leaving the mathematical expressions of the `Gene`s fixed. See `?muthead` for more information.

### Crossover

(To be implemented for system)

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

At the present time, the documentation is the code, unfortunately. But having Julia installed properly, one can execute the files under `tests/` using `julia` to try it out.

## Repository content

* `odes/`: Folder for ODE:s from **[TL]**.
* `tests/`: Folder for testing various aspects of the code.
* `de_package.jl`: File containing the types, functions, and genetic operators.

## TODO/Known issues

Crude TODO-list/feature wish-list

* Implement the rest of the crossover operators listed in **[MF]** and **[E]**
* Also implement their parent-selection methods
* Add support for systems of ODE:s
* Add support for PDE:s (see branch 'partial')
* Add the rest of the ODE:s and NODE:s in **[TL]**
* Document functions and code

Known issues

* It still does not solve all ODE:s (e.g. `ode3`) in reasonable time (if at all).
