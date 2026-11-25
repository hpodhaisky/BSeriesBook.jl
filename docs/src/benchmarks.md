# Benchmarks

The code of BSeriesBook.jl follows the book line by line, and the book's
representation — trees as standard numbers, B-series as plain vectors,
and tables of splittings and partitions that are built once — is also
fast, which matters for research computations at higher orders. This
page compares it with
[BSeries.jl](https://github.com/ranocha/BSeries.jl) 0.1.72 and
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl) 2.27.0
(H. Ranocha, D. I. Ketcheson).

!!! note "Different scope"
    BSeries.jl handles more general objects (arbitrary tree types,
    coloured trees, symbolic coefficients, additive, Rosenbrock and
    multirate methods); BSeriesBook.jl in turn has general linear methods
    and their underlying one-step method. The comparison covers only the
    operations both packages offer.

## Set-up

Script `benchmark/compare.jl` of the repository (BenchmarkTools.jl,
minimum over 5–20 runs). Before each measurement the script checks that
**both packages give the same results**. Machine: Intel Core Ultra 5 225T,
one thread, Julia 1.12.7, BSeriesBook.jl 0.2.0, 25 November 2026.

BSeriesBook times exclude building `Trees(p)` (first row) and the tables
of splittings and partitions, which are built once per `Trees` object on
the first call of `compose` and of `substitute` / `modified_equation`
(their cost is listed below). The factor is the time of BSeries.jl divided
by the time of BSeriesBook.jl (> 1: BSeriesBook.jl is faster).

## Results

### Order 6 (37 trees)

| Task | BSeriesBook.jl | BSeries.jl / RootedTrees.jl | Factor |
|:-----|---:|---:|---:|
| trees + σ + t! | 2.7 µs | 2.3 µs | 0.9 |
| RK B-series, rational | 16.3 µs | 73.1 µs | 4.5 |
| RK B-series, Float64 | 1.5 µs | 13.2 µs | 9 |
| composition, Float64 | 0.6 µs | 50.7 µs | ≈ 85 |
| order of accuracy, rational | 16.4 µs | 72.3 µs | 4.4 |
| modified equation, rational | 359 µs | 561 µs | 1.6 |

### Order 8 (200 trees)

| Task | BSeriesBook.jl | BSeries.jl / RootedTrees.jl | Factor |
|:-----|---:|---:|---:|
| trees + σ + t! | 24.0 µs | 18.1 µs | 0.8 |
| RK B-series, rational | 85.7 µs | 569 µs | 6.6 |
| RK B-series, Float64 | 7.7 µs | 88.9 µs | 12 |
| composition, Float64 | 9.9 µs | 934 µs | 94 |
| order of accuracy, rational | 84.7 µs | 571 µs | 6.7 |
| modified equation, rational | 7.63 ms | 12.9 ms | 1.7 |

### Order 10 (1205 trees)

| Task | BSeriesBook.jl | BSeries.jl / RootedTrees.jl | Factor |
|:-----|---:|---:|---:|
| trees + σ + t! | 599 µs | 156 µs | 0.3 |
| RK B-series, rational | 503 µs | 4.51 ms | 9.0 |
| RK B-series, Float64 | 47.6 µs | 723 µs | 15 |
| composition, Float64 | 257 µs | 18.9 ms | 74 |
| order of accuracy, rational | 501 µs | 4.53 ms | 9.0 |
| modified equation, rational | 348 ms | 768 ms | 2.2 |

### One-off costs and further operations (BSeriesBook.jl)

| | p = 6 | p = 8 | p = 10 |
|:--|---:|---:|---:|
| number of trees | 37 | 200 | 1205 |
| number of splittings (t′, t ∖ t′) | 273 | 3 081 | 37 702 |
| building the table of splittings (once) | 0.02 ms | 0.29 ms | 4.4 ms |
| λ matrix (Algorithm 5), Float64 | 2.1 µs | 38 µs | 1.6 ms |
| `inverse`, Float64 | 0.4 µs | 7.5 µs | 0.21 ms |
| number of partitions (pieces, skeleton) | – | 10 230 | 213 589 |
| building the table of partitions (once) | – | 2.6 ms | 64 ms |
| underlying one-step method, GLM (6.5 h)(ii), rational | – | 4.3 ms | 44 ms |
| `order_of_accuracy`, same GLM (Algorithm 7) | – | 0.12 s | 3.8 s |

BSeries.jl has no general linear methods, so the last row has no
counterpart.

## Where the differences come from

- **Generating trees** is faster in RootedTrees.jl (factor 4 at
  order 10). BSeriesBook.jl also builds the product table `prod`; this is
  the price of the standard numbering, paid once.
- **Elementary weights** are 5–15 times faster: the recursion
  `ψ(t_ℓ ∘ t_r) = ψ(t_ℓ) ⊙ φ(t_r)`, `φ = Aψ`, handles each tree with one
  Hadamard product and one matrix–vector product, and the results live in
  a vector indexed by the standard number instead of a dictionary keyed
  by trees.
- **Composition** is 75–95 times faster: the table of splittings does
  not depend on the two series and is reused, so every composition is a
  single sum over 37 702 entries at order 10. Even including the one-off
  construction of the table (4.4 ms), BSeriesBook.jl is faster at order 10
  than one composition in BSeries.jl.
- **Modified equation** is only about 1.6–2.2 times as fast: both packages spend
  most of the time in `Rational{BigInt}` arithmetic over the partitions,
  not in finding them.
- **Memory** at order 12 (7813 trees): the table of splittings needs about
  43 MB (`Base.summarysize`); a dense λ matrix in Float64 would need
  about 490 MB.

## Reproducing

```
git clone https://github.com/hpodhaisky/BSeriesBook.jl
cd BSeriesBook.jl/benchmark
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
julia --project=. compare.jl
```
