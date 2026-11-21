# BSeriesBook.jl

[![Docs (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://hpodhaisky.github.io/BSeriesBook.jl/stable)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://hpodhaisky.github.io/BSeriesBook.jl/dev)
[![CI](https://github.com/hpodhaisky/BSeriesBook.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/hpodhaisky/BSeriesBook.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](LICENSE)

Julia companion to the book

> J. C. Butcher, *B-Series: Algebraic Analysis of Numerical Methods*,
> Springer Series in Computational Mathematics 55, Springer, 2021
> ([doi:10.1007/978-3-030-70956-3](https://doi.org/10.1007/978-3-030-70956-3));
> 2nd edition in preparation, with H. Podhaisky.

Authors: **John C. Butcher** and **Helmut Podhaisky**.

> **Version 0.1 — initial release, work in progress.** The package
> accompanies a book that is still being revised. It does not yet cover all
> chapters, and names and interfaces may change in later 0.x versions.

The package implements the algorithms of the book **in the book's own
notation**, so that every example can be recomputed: trees are identified
by their standard number (Algorithm 3 of the book), B-series are plain vectors
over this numbering, names follow the book (`σ`, `γ`, `α`, `β`, `λ̂`, `Φ`),
and arithmetic is exact by default.

BSeriesBook.jl is **software for teaching and research**. Beyond the
B-series of Runge–Kutta methods, composition, substitution and modified
equations, it provides the **underlying one-step method of general
linear methods** and the conjugacy test of Algorithm 7 (Chapter 6 of the
book).

For B-series in Julia there is also the package
[BSeries.jl](https://github.com/ranocha/BSeries.jl), together with
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl), by
H. Ranocha and D. I. Ketcheson. It covers further classes of methods
(additive and partitioned Runge–Kutta, Rosenbrock, multirate methods),
coloured trees, symbolic coefficients, symplecticity and energy
preservation. Our test suite uses it as an independent check.

## Installation

```julia
julia> using Pkg; Pkg.add("BSeriesBook")
```

Julia ≥ 1.10. The only run-time dependency is `LinearAlgebra`.

## Example

```julia
julia> using BSeriesBook

julia> T = Trees(5)
Trees: 17 trees up to order 5

julia> [bracket(T, n) for n in 1:4]
4-element Vector{String}:
 "τ"
 "[τ]"
 "[τ²]"
 "[[τ]]"

julia> order_of_accuracy(kutta4(), Trees(8))
4

julia> E = exact_solution(T); half = [E[n] / 2^order(T, n) for n in 1:length(T)];

julia> compose(half, half, T) == E        # two half steps = one step
true
```

A general linear method of Chapter 6 (method (6.5 h)(ii) of the book) and
its underlying one-step method:

```julia
julia> glm = GeneralLinearMethod(
           [0 0 0 0; 3//4 0 0 0; 93//250 -9//125 0 0;
            8881//8640 1069//3240 -1855//5184 0],
           [1 0; 1 1; 1 44//125; 1 691//540],
           [5//54 32//81 250//567 1//14; -19//16 -37//6 1175//336 27//7],
           [1 0; 0 0]);

julia> order_of_accuracy(glm, Trees(7))
5
```

See the [documentation](https://hpodhaisky.github.io/BSeriesBook.jl/stable)
for a tutorial, the table of the standard numbering, the underlying
one-step method, benchmarks against BSeries.jl and the full reference.

## Tests

```julia
julia> using Pkg; Pkg.test("BSeriesBook")
```

The tests compare with the values printed in the book, with RootedTrees.jl
and BSeries.jl, and with mathematical identities.

## Citing

Please cite the book and, if you use the software, the package; see
[`CITATION.bib`](CITATION.bib).

## License

The software is released under the [MIT License](LICENSE),
© 2026 John C. Butcher and Helmut Podhaisky. The license covers the
software only; the text of the book is © Springer Nature.

## Acknowledgements

The conventions for composition and substitution and the test oracles
come from D. I. Ketcheson and H. Ranocha, *Computing with B-series*, ACM
TOMS 49 (2023), [arXiv:2111.11680](https://arxiv.org/abs/2111.11680), and
the packages BSeries.jl and RootedTrees.jl; the substitution law from
P. Chartier, E. Hairer and G. Vilmart, *Algebraic structures of B-series*,
Found. Comput. Math. 10 (2010).
