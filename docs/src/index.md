# BSeriesBook.jl

*Julia companion to J. C. Butcher, **B-Series: Algebraic Analysis of
Numerical Methods**, Springer Series in Computational Mathematics 55,
Springer, 2021 (2nd edition in preparation, with H. Podhaisky).*

!!! warning "Version 0.1 — initial release, work in progress"
    This is the first public version (0.1.0) of the package. It accompanies
    a book that is still being revised, it does not yet cover all chapters,
    and names and interfaces may change in later 0.x versions. Comments and
    bug reports are very welcome on the
    [issue tracker](https://github.com/hpodhaisky/BSeriesBook.jl/issues).

BSeriesBook.jl implements the algorithms of the book **in the book's own
notation**, so that every example can be recomputed:

- trees are identified by their *standard number* `n` from Algorithm 3
  of the book; every tree except τ is the β-product
  `tₙ = t_{L[n]} ∘ t_{R[n]}` (see [Standard numbering](@ref));
- a B-series is a plain vector `a` with `a[n] = a(tₙ)`;
- names follow the book (`σ`, `γ`, `α`, `β`, `λ̂`, `Φ`), with ASCII aliases;
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
[bracket(T, n) for n in 1:8]
```

```@example index
order_of_accuracy(kutta4(), Trees(8))
```

The [Tutorial](@ref) walks through the main functions, the page
[Underlying one-step method](@ref glm) shows the general linear methods
of Chapter 6, and [Benchmarks](@ref) compares the run time with BSeries.jl.

## What is implemented

| Book | Julia |
|:-----|:------|
| Section 2.4, Algorithm 1: number of trees of each order | [`count_trees`](@ref) |
| Section 2.4, Algorithm 2: unrooted and non-superfluous unrooted trees | [`count_unrooted_trees`](@ref) |
| Section 2.5, Algorithm 3: generate trees (standard numbering, `L`, `R`, `prod`) | [`generate_trees`](@ref) |
| Section 2.5, Algorithm 4: generate tree functions (`Bminus`, `factorial`, `symmetry`) | [`tree_functions`](@ref) |
| Section 2.5: σ, t!, α, β; β-product; B⁺ form | [`Trees`](@ref), [`σ`](@ref), [`γ`](@ref), [`α`](@ref), [`β`](@ref), [`beta_product`](@ref), [`children`](@ref), [`bracket`](@ref) |
| Section 3.5: exact solution E(t) = 1/t! | [`exact_solution`](@ref) |
| Section 3.9, Algorithms 5 and 6: evaluate λ (and its efficient form λ̂) | [`lambda_matrix`](@ref), [`lambda_hat`](@ref) |
| Section 3.9: composition and inverse | [`compose`](@ref), [`inverse`](@ref) |
| Section 3.9: substitution (Chartier–Hairer–Vilmart) | [`substitute`](@ref) |
| Backward error analysis | [`modified_equation`](@ref), [`modifying_integrator`](@ref) |
| Section 3.6: elementary weights, order of RK methods | [`elementary_weights`](@ref), [`order_of_accuracy`](@ref) |
| Section 5.3: stability function | [`stability_polynomial_coefficients`](@ref) |
| Section 6.6, Algorithm 7: are two central B-series conjugate? | [`conjugacy_step`](@ref), [`conjugate_order`](@ref) |
| Section 6.4: general linear methods | [`GeneralLinearMethod`](@ref), [`stability_matrix`](@ref) |
| Section 6.6: partitioned form (6.6 a), underlying one-step method (6.6 d), order | [`partitioned_form`](@ref), [`underlying_one_step`](@ref) |

Not yet implemented (planned for later versions): elementary differentials
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
  note   = {Version 0.1},
  url    = {https://github.com/hpodhaisky/BSeriesBook.jl}
}
```

## License

MIT License, © 2026 John C. Butcher and Helmut Podhaisky. The license
covers the software only; the text of the book is © Springer Nature.
