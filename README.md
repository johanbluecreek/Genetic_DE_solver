# Genetic DE solver

This is a differential equation solver that finds analytical solutions using genetic programming. This is written in [Julia](https://github.com/JuliaLang/julia) and relies on the [Calculus](https://github.com/johnmyleswhite/Calculus.jl) package for symbolic derivatives, and fitness is calculated with the help of [meval-rs](https://github.com/rekka/meval-rs), a [Rust](https://www.rust-lang.org) library.

## Introduction

Finding analytical solutions to differential equations (DEs) is quite a hard problem. There are some solvers out there, see e.g. [SymPy/solvers](http://docs.sympy.org/latest/tutorial/solvers.html#solving-differential-equations), that are able to give you solutions to particular DEs if the form of them are simple enough. A general DE however, even if it does have a solution in terms of elementary functions, cannot be solved by such solvers.

This is where the present solver would be relevant. Say your choice of solver does not provide you with a solution, perhaps there still is one. By generating random starting candidate solutions, this solver applies genetic operations -- mutation and crossover -- to attempt to generate a full solution.

The more correct denotation of this project would be *frame-work* and not *solver*. There are a multitude of parameters that can be given as input; from the set of functions the algorithm will use to the various probabilities for mutation or crossover to happen. Hence there are a multitude of solvers that can be implemented using this frame-work. The aim here is to both provide the frame-work and some pre-defined solvers that can be used directly.

To get started, take a look at the *Usage* section. A description of the frame-work is outlined under the *Project overview* section, and a description of features and goals can be found in the *Status* section.

The inspiration for this project is taken from

* **[TL]** TSOULOS, Ioannis G. et LAGARIS, Isaac E. "*Solving differential equations with genetic programming*". Genetic Programming and Evolvable Machines, 2006, vol. 7, no 1, p. 33-54.

the following books have been used as a reference

* **[MF]** "*How to Solve It: Modern Heuristics*" by Zbigniew Michalewicz, David
B. Fogel
* **[E]** "*Computational Intelligence: An Introduction*" by Andries P.
Engelbrecht

other references

* **[TGG]** TSOULOS, Ioannis G. et GAVRILIS, Dimitris et GLAVAS, Euripidis "*Solving differential equations with constructed neural networks*". Neurocomputing 72 (2009) 2385-2391.

## Usage

Usage documentation is very lacking at the moment. But, it has been tried on a machine with the following software

 * Julia v0.6.1
 * Julia packages Calculus and Iterators installed (`Pkg.add()`)
 * Rust v1.19.0
 * Cargo v0.21.0

The code is mostly Julia, however to circumvent some bugs in Julia there is a Rust library present that performs the fitness calculation. Before attempting to run any of the code, you should build that library

```$ cargo build```

from the root of this repository.

To start solving differential equations, there is yet not much user-friendly instructions nor documentation, but you can take a look at `src/solvers/TL_solver.jl` to begin with, which is the most mature solver provided in this frame-work as of now. Just try to run it, e.g.

```$ julia src/solvers/TL_solver.jl```

### Suggested work-flow

These are still early days for this project, and if you want to try it out, I would suggest to proceed in the following manner:

 * Create a file containing the (system of) differential equation(s) you want to solve. See the files under `des/` for examples of how to write such files.
 * Create a second file for your solver. Load the necessary functions by including `src/GeneticDESolver.jl`, and include your file containing the differential equation. See the files under `solvers/` for examples.

## Project overview

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

The above table depicts the structure of a solution to a system of `n` functions, where each function is made up by `m` sub-expressions.

It should be noted that the number of differential equations is independent from the number of unknown functions `n`. You can provide the algorithm with more (*over-determined*) or fewer (*under-determined*) differential equations (or even algebraic constraints) than functions.

The sub-expression number `m` is independent of any details of the problem, and can be chosen at will. It should however be noted that larger `m` will have negative effect on the time it takes to execute each iteration. See **[MF]** regarding the introduction of sub-expressions.

**[TL]** Uses a [BNF-grammar](https://en.wikipedia.org/wiki/Backus%E2%80%93Naur_form) while the current project has uses [Polish notation](https://en.wikipedia.org/wiki/Polish_notation) for representing expressions.

### Mutation

There are several mutation operators implemented, and they are divided in two parts.

**Gene mutation:** The main form of mutations acts on the `Gene`s and changes the mathematical expressions they give rise to. See `?mutate` for further explanation on this, and the Status section below.

**Chromosome mutation:** A supplementary form of mutation that changes how the `Gene`s are combined to form a `Chromosome`, leaving the mathematical expressions of the `Gene`s fixed. See `?muthead` for more information, and the Status section below.

### Crossover

There are some crossover operators implemented. These are interactions between `Individual`s such that they produce new solutions as offspring. See `?crossover` for more information, and the Status section below.

### Fitness calculation

The fitness if calculated similar to that of **[TL]**. The fitness is here divided in three parts

* **error:** Sum of squares of the differential equations evaluated at a set number of uniformly distributed points over a given interval. (Also in **[TL]**)
* **penalty:** Absolute value of the mismatch for the boundary condition. (Also in **[TL]**)
* **shape:** This is similar to error, but evaluated for derivatives of the differential equation, with a decay-factor for each derivative. (Not in **[TL]**)

both penalty and shape can be turned off at your preference. The fitness is then simply given as the sum of these parts.

## Status

Here is an incomplete list of things that are/should be implemented.

### Large scale stuff

  * Meta-programming overhaul. ✗

    Presently the code does not make use of Julia's built in meta-programming features very much. This should be improved.

  * Expand on the function basis. ✗

    The function set is currently bound by `meval-rs`, which supports only functions in the Rust std library. So first thing would be...

    - Extend functions supported by `meval-rs`.

    Important functions would be, for example, [Generalised Hypergeometric Functions](https://en.wikipedia.org/wiki/Generalized_hypergeometric_function).

  * Complement with neural-network option, as in **[TGG]**. ✗

### Capabilities

Overall capabilities

|         | ODE | SODE | PDE | SPDE |
|---------|:---:|:----:|:---:|:----:|
| Works   |  ✓  |  ✓   |  ✓  |  ?   |
| Problem |     |      |     |  not tested    |

### Fitness and similar

|         | ODE | SODE | PDE | SPDE |
|---------|:---:|:----:|:---:|:----:|
| Error   |  ✓  |  ✓   |  ✓  |  ?   |
| Penalty |  ✓  |  ✓   |  ✓  |  ?   |
| Shape   |  ✗  |  ✗   |  ✗  |  ✗   |

✓ = implemented, ✗ = not implemented, ? = not tested and/or not intensionally implemented

### Genetic operators

Table of implemented methods for the genetic operators

| `mutate` | `crossover` | `muthead`| `p_select`       |
| -------- | ----------- | -------- | ---------------- |
| change   | safe1point  | scramble | tournament       |
| swap     | 1point      | jump     | random           |
| grow     | 2point      | combo    |
| trunc    | random      |
| random   |

You can find documentation on each function and references to documentation of each method by the name presented above on the first row.

### Provided differential equations

From **[TL]** the following differential equations are provided with this repository

| ODE | NLODE | SODE | PDE  |
| --- | ----- | ---- | ---- |
| 1   |  3    |  1   | 1    |
| 3   |       |      |      |

Other differential equations (no boundary conditions)

| ODE | NLODE | SODE | PDE       | SPDE |
| --- | ----- | ---- | --------- | ---- |
|||| [Liouville](https://en.wikipedia.org/wiki/Liouville%27s_equation) ||


## Contributing

Any contribution would be appreciated, you can start off small, e.g.

* Write the files for any of the differential equation in **[TL]**.

try it out and tune it, e.g.

* Find out optimal settings for mutation (rate, selection, sizes, operations), crossover, etc. for the different differential equations, for example `ode3.jl`.

or help developing, for example

* See the TODO below
* Read TODO, XXX, and FIXME comments in the code

## Repository content

* `src/`: Source files.
  - `src/GeneticDESolver.jl`: Main file that makes all functions available. Type `include("src/GeneticDESolver.jl")` to load.
  - `src/lib.rs`: Rust library to assist fitness calculation.
  - `src/genops/`: Genetic operators.
* `des/`: Folder containing pre-defined differential equations (and systems of), primarily for benchmarking purposes and tests.
* `tests/`: Folder containing tests of various aspects of the code (OLD)
* `solvers/`: Folder containing pre-defined solvers. (TODO -- move `TL_solver.jl` here)

## TODO

Crude TODO-list/feature wish-list

* Implement the rest of the crossover operators listed in **[MF]** and **[E]**
* Also implement their parent-selection methods
* Add the rest of the differential equations in **[TL]**
* Complete the items marked with "✗" of "?" in the Status section.
* Make a pedagogical introduction in a jupyter-notebook for how all this works.
* Construct benchmarks
* Construct a proper test suite.
