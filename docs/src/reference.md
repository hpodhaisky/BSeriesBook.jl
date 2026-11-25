# Reference

The functions are grouped by the chapters of the book.

```@meta
CurrentModule = BSeriesBook
```

```@docs
BSeriesBook
```

## Data types

| Type | What it is | Book |
|:-----|:-----------|:-----|
| [`Trees`](@ref) | all rooted trees up to a given order, with their tree functions; `Trees(p_max)` | 2.5 |
| [`Tree`](@ref) | one rooted tree t_n: its standard number `n` together with its `Trees`; `T[n]` | 2.5 |
| [`Forest`](@ref) | a multiset of trees, e.g. B⁻(t) = `children(t)` | 2.1 |
| [`Splitting`](@ref) | a pair (t′, t ∖ t′): subtree containing the root, forest of removed parts; `splittings(t)` | 2.8 |
| [`Partition`](@ref) | a partition of t into pieces by cutting edges, with its skeleton; `tree_partitions(t)` | 3.9 |
| [`AbstractBSeries`](@ref) | supertype of the three B-series types; a vector of coefficients a(t), indexed by trees | 3.4 |
| [`GeneralBSeries`](@ref) | a B-series in 𝖡* with an arbitrary value a(∅) (stored); a vector space | 3.4 |
| [`CentralBSeries`](@ref) | a B-series in 𝖡 with a(∅) = 1, e.g. of a one-step method; the B-group under composition | 3.4 |
| [`VectorFieldBSeries`](@ref) | a B-series in 𝖡⁰ with a(∅) = 0; h⁻¹ B_h(a) is a vector field, e.g. a modified equation | 3.4 |
| [`RungeKuttaMethod`](@ref) | a Runge–Kutta tableau (A, b, c) | 1, 3.6 |
| [`GeneralLinearMethod`](@ref) | a general linear method (A, U, B, V) | 6.4 |

[`empty_tree`](@ref) (alias `∅`) is the index of
a(∅) in a B-series.

## Chapter 2: trees

```@docs
Trees
Tree
number
order
symmetry
factorial(::Tree)
alpha
beta
trees_of_order
beta_product
Forest
children
bplus
bracket
Splitting
splittings
count_trees
```

## Chapter 3: B-series

```@docs
AbstractBSeries
GeneralBSeries
CentralBSeries
VectorFieldBSeries
empty_tree
coefficients
Base.zero(::AbstractBSeries)
exact_solution
identity_series
rescale
slope
lambda_matrix
lambda_hat
compose
inverse
```

## Chapter 3: substitution and modified equations

```@docs
Partition
tree_partitions
substitute
modified_equation
modifying_integrator
```

## Chapters 1, 3, 5: Runge–Kutta methods

```@docs
RungeKuttaMethod
stages
elementary_weights
bseries
order_of_accuracy
stability_polynomial_coefficients
euler
runge_midpoint
kutta4
kutta38
gauss2
butcher6a
```

## Chapter 6: general linear methods

```@docs
conjugacy_step
conjugate_order
GeneralLinearMethod
inputs
is_partitioned_form
partitioned_form
underlying_one_step
```

## Internals

Not exported; use them with the module prefix, e.g.
`BSeriesBook.generate_trees(4)`.  The building blocks of [`Trees`](@ref)
are the tables of Algorithms 3 and 4 of the book.

```@docs
generate_trees
TreeTable
ProductTable
tree_functions
TreeFunctions
level_sequence
join_roots
count_unrooted_trees
stage_weights
stability_matrix
```

## Index

```@index
```
