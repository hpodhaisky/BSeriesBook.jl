# BSeriesBook.jl

*Julia companion to J. C. Butcher, **B-Series: Algebraic Analysis of
Numerical Methods**, Springer Series in Computational Mathematics 55,
Springer, 2021 (2nd edition in preparation, with H. Podhaisky).*

!!! warning "Version 0.2 — work in progress"
    The package accompanies a book that is still being revised, it does
    not yet cover all chapters, and names and interfaces may change in
    later 0.x versions (see `NEWS.md` in the repository for the changes
    from 0.1). Comments and
    bug reports are very welcome on the
    [issue tracker](https://github.com/hpodhaisky/BSeriesBook.jl/issues).

BSeriesBook.jl implements the algorithms of the book **in the book's own
notation**, so that every example can be recomputed:

- trees are identified by their *standard number* `n` from Algorithm 3
  of the book; every tree except τ is the β-product
  `tₙ = t_{L[n]} ∘ t_{R[n]}` (see the [Tutorial](@ref standard-numbering));
- a B-series is a vector `a` indexed by trees, `a[t] = a(t)`; its type records
  a(∅) as in the book: [`GeneralBSeries`](@ref) (𝖡\*),
  [`CentralBSeries`](@ref) (𝖡, a(∅) = 1) and
  [`VectorFieldBSeries`](@ref) (𝖡⁰, a(∅) = 0);
- names follow the book, in ASCII: `symmetry` (σ), `factorial` (t!),
  `alpha`, `beta`, `lambda_hat` (λ̂), `slope`; the Greek symbols `σ`,
  `α`, `β` and `∅` are also available as aliases;
- arithmetic is exact (`Rational{BigInt}`) by default;
- algorithms, sections, equations and tables are cited with their numbers
  in the printed book (1st edition), e.g. "Algorithm 3", "(6.6 d)".

BSeriesBook.jl is **software for teaching and research**. Besides the
B-series of Runge–Kutta methods, composition, substitution and modified
equations, it contains material of the book that other packages do not
cover, most importantly the
[underlying one-step method of a general linear method](@ref glm)
and the conjugacy test of Algorithm 7 (Chapter 6). The data structures of the
book are also fast; see [Benchmarks](@ref).

For B-series in Julia there is also the package
[BSeries.jl](https://github.com/ranocha/BSeries.jl), together with
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl), by
H. Ranocha and D. I. Ketcheson (see Ketcheson and Ranocha, *Computing
with B-series*, ACM TOMS 2023). It covers further classes of methods
(additive and partitioned Runge–Kutta, Rosenbrock, multirate methods),
coloured trees, symbolic coefficients, symplecticity and energy
preservation. BSeriesBook.jl is tested against it as an independent
reference. Further software for B-series is listed under
[Literature](@ref).

## Installation

BSeriesBook.jl is registered in the Julia
[General registry](https://github.com/JuliaRegistries/General),
so it is installed with the package manager:

```julia
julia> using Pkg

julia> Pkg.add("BSeriesBook")
```

Julia 1.10 or newer is required. The only run-time dependency is the
standard library `LinearAlgebra`.

## First steps

```@example index
using BSeriesBook
T = Trees(5)
```

```@example index
[bracket(t) for t in T if order(t) <= 4]
```

```@example index
order_of_accuracy(kutta4(), Trees(8))
```

The [Tutorial](@ref) explains the conventions and walks through the main
computations, chapter by chapter, up to the
[underlying one-step method](@ref glm) of Chapter 6. The
[Reference](@ref) lists all functions in the order of the book;
[Benchmarks](@ref) compares the run time with BSeries.jl.

## Not yet implemented

Planned for later versions: elementary differentials
and numerical evaluation of B-series, symplecticity and energy
preservation, G-symplectic methods of Chapter 7, the integration methods
and the B-group of Chapter 4.

## Citing

If you use the package, please cite the book and the package:

```bibtex
@book{Butcher2021,
  author    = {Butcher, John C.},
  title     = {B-Series: Algebraic Analysis of Numerical Methods},
  series    = {Springer Series in Computational Mathematics},
  volume    = {55},
  publisher = {Springer},
  address   = {Cham},
  year      = {2021},
  doi       = {10.1007/978-3-030-70956-3}
}

@misc{BSeriesBook,
  author = {Butcher, John C. and Podhaisky, Helmut},
  title  = {{BSeriesBook.jl}: {J}ulia companion to the book
            ``{B}-Series: Algebraic Analysis of Numerical Methods''},
  year   = {2026},
  note   = {Version 0.2},
  url    = {https://github.com/hpodhaisky/BSeriesBook.jl}
}
```

## License

MIT License, © 2026 John C. Butcher and Helmut Podhaisky. The license
covers the software only; the text of the book is © Springer Nature.
