# Conventions

## Trees

Every tree except `τ` is the β-product of a unique pair,
`tₙ = t_{L[n]} ∘ t_{R[n]}` (graft `t_{R[n]}` onto the root of `t_{L[n]}`),
with `R[n] ≤ R[L[n]]`. Algorithm 3 of the book ([`generate_trees`](@ref))
builds the tables `order`, `first`, `last`, `L`, `R` and `prod` with the
book's 1-based indices.

In the book a tree *is* its standard number. In the package a tree is a
[`Tree`](@ref): `t = T[n]` carries the number `n` = [`number`](@ref)`(t)`
and its `Trees`, but it is not an integer, so a tree cannot be confused
with an order, a count or an index. Functions on trees take a `Tree`
(`order(t)`, `factorial(t)`), and a B-series is indexed by trees,
`a[t]`. The type costs nothing at run time: a loop over `Tree`s compiles
to the same code as a loop over the numbers. The tables (`L`, `R`,
`prod`, [`splittings`](@ref), [`tree_partitions`](@ref)) stay tables of
standard numbers, as in the book's algorithms.

| Representation | Function |
|:---------------|:---------|
| standard number `n` | `T[n]`, [`number`](@ref)`(t)`; [`trees_of_order`](@ref)`(T, p)`, `for t in T` |
| β-product `s ∘ t` | [`beta_product`](@ref)`(s, t)` |
| B⁺ form `[t₁ ⋯ t_m]` | [`children`](@ref), inverse [`bplus`](@ref), text [`bracket`](@ref) |
| level sequence | [`level_sequence`](@ref), for conversion to RootedTrees.jl |

The product table `prod` is stored sparsely: only pairs with
`order[ℓ] + order[r] ≤ p_max` exist (a dense table would have about 61
million entries at order 12).

## B-series

A B-series is a mapping a: T# → ℝ, T# = {∅} ∪ T, and stands for the
series

```math
B_h(a)\,y = a(\emptyset)\,y + \sum_{t} \frac{h^{|t|}\,a(t)}{\sigma(t)}\,F(t)(y),
```

the same normalisation as in BSeries.jl; the exact solution is
`E(t) = 1/t!`. In the package it is a vector over the standard numbering,
`a[n] = a(tₙ)`, truncated at the order of its [`Trees`](@ref); the value
on the empty tree is `a[empty_tree]` (or `a[∅]`). The type records a(∅), following the three
sets of Section 3.4 of the book:

| Book | a(∅) | Type | Structure |
|:-----|:-----|:-----|:----------|
| 𝖡\* | stored | [`GeneralBSeries`](@ref) | vector space |
| 𝖡 | 1 (not stored) | [`CentralBSeries`](@ref) | group under composition (the B-group) |
| 𝖡⁰ | 0 (not stored) | [`VectorFieldBSeries`](@ref) | linear space; `h⁻¹ B_h(a)` is a vector field |

The type is decided by the structure of a computation, never by the
value of a(∅): `a - b` of two central series is in 𝖡⁰, `a + v` with
`v` in 𝖡⁰ is central, `a + b` and `2a` are general. A conversion such
as `CentralBSeries(g)` checks the value of a(∅). Equality `==` and `≈`
compare a(∅) as well; indexing with a range, broadcasting and other
generic vector operations act on the coefficients a(tₙ) only and return
plain vectors ([`coefficients`](@ref) gives the stored vector itself).

## Composition

[`compose`](@ref)`(a, b)`, also written `a * b`, is the book's product
`ab`: **first `a`, then `b`** (Theorem 3.9), a map 𝖡 × 𝖡\* → 𝖡\*,

```math
(ab)(\emptyset) = b(\emptyset), \qquad
(ab)(t) = b(\emptyset)\,a(t) + \sum_{\emptyset < t' \le t} a(t \setminus t')\, b(t'),
```

where `t'` runs over the subtrees of `t` that contain the root, and
`a(F)` is the product of `a` over the trees of the forest `F`. This
agrees with `BSeries.compose(a, b)`. The first factor must be central
(otherwise `B_h(b) ∘ B_h(a)` is not a B-series); the result has the type
of `b`. In 𝖡 the inverse is [`inverse`](@ref)`(a)` = `inv(a)`, the
identity `one(a)` = [`identity_series`](@ref)`(T)`, and `a^n` is defined
for all integers `n`.

The sum is taken from the table of [`splittings`](@ref), built once per
`Trees` object by the β-product recursion and reused for every
composition.

## The map slope

The B-series of `hf(y₀)` is D ∈ 𝖡⁰ (D(τ) = 1, D(t) = 0 otherwise), and
for a central `a` the series of `hf(B_h(a)y₀)` is `aD`, with
`(aD)([t₁ ⋯ t_m]) = a(t₁) ⋯ a(t_m)`. [`slope`](@ref)`(y)` computes
`yD`; `slope(y, aii)` solves the implicit relation of a diagonally
implicit stage, returning `dw = wD` with `w = y + aii·dw`, in one pass
over the standard numbering. `slope(identity_series(T))` is D.

## Substitution

[`substitute`](@ref)`(b, a)` is `b ⋆ a` for `b` in 𝖡⁰ (argument order as
in `BSeries.substitute`): the method `a` applied to the vector field
`h⁻¹ B_h(b)`,

```math
(b \star a)(t) = \sum_{\text{partitions of } t} b(\text{piece}_1)\cdots b(\text{piece}_m)\; a(\text{skeleton}),
```

with `(b ⋆ a)(∅) = a(∅)`; the result has the type of `a`. The
partitions with their skeletons ([`tree_partitions`](@ref)) are built by
the same β-product recursion as the splittings.

## Number types

Functions are generic in the number type. The tableaux of the book are
given as rationals; [`exact_solution`](@ref) returns `Rational{BigInt}`,
and [`GeneralLinearMethod`](@ref) widens integer and rational
coefficients to `Rational{BigInt}`, so that no overflow occurs at high
orders. Floating-point tableaux stay floating point (use the `atol`
keyword of [`order_of_accuracy`](@ref)).
