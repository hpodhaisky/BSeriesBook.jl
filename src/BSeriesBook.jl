"""
    BSeriesBook

Julia companion to J. C. Butcher, *B-Series: Algebraic Analysis of
Numerical Methods*, Springer Series in Computational Mathematics 55
(2nd edition in preparation).

The package follows the book's notation and algorithms: trees are
identified by their standard number (Algorithm 3), and B-series are
vectors over this numbering whose type records a(∅): `GeneralBSeries`
(𝖡*), `CentralBSeries` (𝖡, a(∅) = 1) and `VectorFieldBSeries`
(𝖡⁰, a(∅) = 0).  It is software for teaching and research;
in addition to the B-series of Runge–Kutta methods it provides the
underlying one-step method of general linear methods (Chapter 6).

For B-series in Julia see also BSeries.jl and RootedTrees.jl
(H. Ranocha, D. I. Ketcheson), which the test suite uses as an
independent reference.
"""
module BSeriesBook

include("trees.jl")
include("enumeration.jl")
include("bseries.jl")
include("runge_kutta.jl")
include("substitution.jl")
include("glm.jl")
include("underlying.jl")

# Not exported (use with the prefix `BSeriesBook.`): the building blocks
# of `Trees` (`generate_trees`, `tree_functions`, `TreeTable`,
# `TreeFunctions`, `ProductTable`), `level_sequence`, `join_roots`,
# `count_unrooted_trees`, `stage_weights`, `stability_matrix`.

# Chapter 2
export Trees, Tree, number
export order, symmetry, alpha, beta, σ, α, β   # and Base.factorial(t)
export trees_of_order, beta_product, children, bplus, Forest
export Splitting, splittings
export bracket
export count_trees

# Chapter 3
export AbstractBSeries, GeneralBSeries, CentralBSeries, VectorFieldBSeries
export empty_tree, ∅, coefficients
export exact_solution, identity_series, rescale
export lambda_matrix, lambda_hat, compose, inverse, slope
export Partition, tree_partitions
export substitute, modified_equation, modifying_integrator

# Chapters 1, 3, 5
export RungeKuttaMethod, stages, elementary_weights, bseries
export order_of_accuracy, stability_polynomial_coefficients
export euler, runge_midpoint, kutta4, kutta38, gauss2
export butcher6a

# Chapter 6
export conjugacy_step
export GeneralLinearMethod, inputs
export is_partitioned_form, partitioned_form, underlying_one_step
export conjugate_order

end # module
