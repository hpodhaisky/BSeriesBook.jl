# News

## Version 0.2.0

Breaking release: trees and B-series have types.

- A tree is a `Tree`, no longer an integer: `T[n]` is the tree with
  standard number `n`, `number(t)` gives `n`, and iterating over `T`
  gives all trees. Tree functions take the tree: `order(t)`,
  `symmetry(t)`, `factorial(t)`, `alpha(t)`, `beta(t)`, `bracket(t)`,
  `level_sequence(t)`, `children(t)`, `beta_product(s, t)`,
  `join_roots(s, t)`; `trees_of_order(T, p)` returns trees;
  `tree_number(T, kids)` is now `bplus(T, kids)`. B-series are indexed by
  trees, `a[t]`, and the function constructors receive trees,
  `CentralBSeries(t -> 1 // factorial(t), T)`. No run-time cost: a loop
  over trees compiles to the same code as a loop over numbers. The
  tables (`L`, `R`, `prod`, splittings, partitions) stay tables of
  standard numbers.

- B-series are no longer plain vectors but have one of three types, the
  sets of Section 3.4 of the book: `GeneralBSeries` (𝖡\*, a(∅) stored),
  `CentralBSeries` (𝖡, a(∅) = 1) and `VectorFieldBSeries` (𝖡⁰, a(∅) = 0),
  all subtypes of `AbstractBSeries <: AbstractVector`. Each stores its
  `Trees`, so the argument `T` is no longer passed to functions of
  B-series. `a[∅]` (or `a[empty_tree]`) is a(∅); `coefficients(a)` is the
  vector of the a(tₙ).
- `+`, `-`, scalar multiplication and division return the type that
  follows from a(∅): the difference of two central series lies in 𝖡⁰.
  `==` and `≈` compare a(∅) as well.
- `compose(a, b)` (also `a * b`) requires a central first factor and
  returns the type of `b`; the keyword `b0` is gone. `inverse(a)`
  (also `inv(a)`), `one(a)` and `a^n` are defined for central series.
- `substitute(b, a)` requires `b` in 𝖡⁰ and returns the type of `a`.
  `modified_equation(a)` and `modifying_integrator(a)` return a
  `VectorFieldBSeries`.
- `exact_solution`, `identity_series` and `bseries(rk, T)` return a
  `CentralBSeries`; `elementary_weights` still returns a plain vector.
- `lambda_matrix(a, p_top)`, `lambda_hat(a, p_top)`: `p_top` defaults to
  the order of the trees of `a`.
- `conjugacy_step(a, b, p)` returns `(; x, conjugated)` on success and `nothing`
  otherwise (instead of `(true, a′, x)` / `(false,)`);
  `conjugate_order(a, b)`.
- `underlying_one_step(glm, T)` returns `E` as a `CentralBSeries` and the
  starting method as `S`, a vector of `GeneralBSeries` (one per input,
  `S[i][empty_tree] = u[i]`), instead of the vectors `ζ`; the fields `ξ`, `η`,
  `ηD` are now `xi`, `eta`, `etaD`.
- New: `slope(y, aii = 0)`, the series `yD` of `hf(B_h(y)y₀)`, and for
  `aii ≠ 0` the solution of the implicit relation `w = y + aii·wD` of a
  diagonally implicit stage (the function f of J. C. Butcher).
- New: `CentralBSeries{S}(a)` (and likewise for the other types) changes
  the number type, e.g. `CentralBSeries{BigFloat}(exact_solution(T))`.
- New: `rescale(a, θ)`, the B-series for the step size θh
  (`rescale(E, θ)` = E^(θ)); `order_of_accuracy(a::CentralBSeries)`.
- B-series are displayed with the standard number and the bracket form
  of each tree.
- ASCII names first: `symmetry`, `alpha`, `beta`, `empty_tree`; the
  Greek symbols `σ`, `α`, `β`, `∅` remain as aliases. The tree factorial
  t! is `factorial(t)`, a method of `Base.factorial` (the book's name;
  `γ` and `density` are gone).

## Version 0.1.0

Initial release.
