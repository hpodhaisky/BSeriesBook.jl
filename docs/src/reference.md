# Reference

The functions are grouped by the chapters of the book.

```@meta
CurrentModule = BSeriesBook
```

```@docs
BSeriesBook
```

## Chapter 2: trees

```@docs
Trees
Tree
number
generate_trees
TreeTable
ProductTable
tree_functions
TreeFunctions
order
symmetry
factorial(::Tree)
alpha
beta
trees_of_order
beta_product
children
bplus
level_sequence
bracket
Splitting
splittings
count_trees
count_unrooted_trees
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
join_roots
substitute
modified_equation
modifying_integrator
```

## Chapters 1, 3, 5: Runge–Kutta methods

```@docs
RungeKuttaMethod
stages
stage_weights
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
stability_matrix
is_partitioned_form
partitioned_form
underlying_one_step
```

## Index

```@index
```
