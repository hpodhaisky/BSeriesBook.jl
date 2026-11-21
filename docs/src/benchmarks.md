# Benchmarks

The code of BSeriesBook.jl follows the book line by line, and the book's
representation — trees as standard numbers, B-series as plain vectors,
and tables of splittings and partitions that are built once — is also
fast, which matters for research computations at higher orders. This
page compares it with
[BSeries.jl](https://github.com/ranocha/BSeries.jl) 0.1.72 and
[RootedTrees.jl](https://github.com/SciML/RootedTrees.jl) 2.26.0
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
one thread, Julia 1.12, 18 November 2026.

BSeriesBook times exclude building `Trees(p)` (first row) and the tables
of splittings and partitions, which are built once per `Trees` object on
the first call of `compose` and of `substitute` / `modified_equation`
(their cost is listed below). The factor is the time of BSeries.jl divided
by the time of BSeriesBook.jl (> 1: BSeriesBook.jl is faster).

## Results

### Order 6 (37 trees)

| Task | BSeriesBook.jl | BSeries.jl / RootedTrees.jl | Factor |
|:-----|---:|---:|---:|
| trees + σ + t! | 2.9 µs | 2.4 µs | 0.8 |
| RK B-series, rational | 17.0 µs | 76.5 µs | 4.5 |
| RK B-series, Float64 | 1.7 µs | 12.6 µs | 7 |
| composition, Float64 | 0.6 µs | 54.5 µs | ≈ 90 |
| order of accuracy, rational | 17.0 µs | 76.5 µs | 4.5 |
| modified equation, rational | 353 µs | 536 µs | 1.5 |

### Order 8 (200 trees)

| Task | BSeriesBook.jl | BSeries.jl / RootedTrees.jl | Factor |
|:-----|---:|---:|---:|
| trees + σ + t! | 29.7 µs | 17.4 µs | 0.6 |
| RK B-series, rational | 86.1 µs | 575 µs | 6.7 |
| RK B-series, Float64 | 8.2 µs | 94.1 µs | 11 |
| composition, Float64 | 9.3 µs | 1.06 ms | 114 |
| order of accuracy, rational | 86.0 µs | 572 µs | 6.6 |
| modified equation, rational | 7.41 ms | 13.0 ms | 1.8 |

### Order 10 (1205 trees)

| Task | BSeriesBook.jl | BSeries.jl / RootedTrees.jl | Factor |
|:-----|---:|---:|---:|
| trees + σ + t! | 567 µs | 157 µs | 0.3 |
| RK B-series, rational | 507 µs | 4.52 ms | 8.9 |
| RK B-series, Float64 | 47.8 µs | 718 µs | 15 |
| composition, Float64 | 236 µs | 19.4 ms | 82 |
| order of accuracy, rational | 507 µs | 4.57 ms | 9.0 |
| modified equation, rational | 381 ms | 735 ms | 1.9 |

### One-off costs and further operations (BSeriesBook.jl)

| | p = 6 | p = 8 | p = 10 |
|:--|---:|---:|---:|
| number of trees | 37 | 200 | 1205 |
| number of splittings (t′, t ∖ t′) | 273 | 3 081 | 37 702 |
| building the table of splittings (once) | 0.06 ms | 0.75 ms | 9.4 ms |
| λ matrix (Algorithm 5), Float64 | 3.3 µs | 75 µs | 1.7 ms |
| `inverse`, Float64 | 0.5 µs | 8.8 µs | 0.33 ms |
| number of partitions (pieces, skeleton) | – | 10 230 | 213 589 |
| building the table of partitions (once) | – | 2.6 ms | 0.11 s |
| underlying one-step method, GLM (6.5 h)(ii), rational | – | 4.8 ms | 97 ms |

BSeries.jl has no general linear methods, so the last row has no
counterpart.

## Where the differences come from

- **Generating trees** is faster in RootedTrees.jl (factor 2–4 at
  order 10). BSeriesBook.jl also builds the product table `prod`; this is
  the price of the standard numbering, paid once.
- **Elementary weights** are 5–15 times faster: the recursion
  `ψ(t_ℓ ∘ t_r) = ψ(t_ℓ) ⊙ φ(t_r)`, `φ = Aψ`, handles each tree with one
  Hadamard product and one matrix–vector product, and the results live in
  a vector indexed by the standard number instead of a dictionary keyed
  by trees.
- **Composition** is 80–110 times faster: the table of splittings does
  not depend on the two series and is reused, so every composition is a
  single sum over 37 702 entries at order 10. Even including the one-off
  construction of the table (9.4 ms), BSeriesBook.jl is faster at order 10
  than one composition in BSeries.jl.
- **Modified equation** is only about twice as fast: both packages spend
  most of the time in `Rational{BigInt}` arithmetic over the partitions,
  not in finding them.
- **Memory** at order 12 (7813 trees): the table of splittings needs about
  100 MB; a dense λ matrix would need about 470 MB.

## Reproducing

```
git clone https://github.com/hpodhaisky/BSeriesBook.jl
cd BSeriesBook.jl/benchmark
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
julia --project=. compare.jl
```
