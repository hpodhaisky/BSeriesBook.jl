# Conventions

## Trees

Every tree except `τ` is the β-product of a unique pair,
`tₙ = t_{L[n]} ∘ t_{R[n]}` (graft `t_{R[n]}` onto the root of `t_{L[n]}`),
with `R[n] ≤ R[L[n]]`. Algorithm 3 of the book ([`generate_trees`](@ref))
builds the tables `order`, `first`, `last`, `L`, `R` and `prod` with the
book's 1-based indices; a tree *is* its standard number `n`.

| Representation | Function |
|:---------------|:---------|
| standard number `n` | basic; `trees_of_order(T, p)` is `first[p]:last[p]` |
| β-product `t_ℓ ∘ t_r` | [`beta_product`](@ref)`(T, ℓ, r)` |
| B⁺ form `[t_{k₁} ⋯ t_{k_m}]` | [`children`](@ref), inverse [`tree_number`](@ref), text [`bracket`](@ref) |
| level sequence | [`level_sequence`](@ref), for conversion to RootedTrees.jl |

The product table `prod` is stored sparsely: only pairs with
`order[ℓ] + order[r] ≤ p_max` exist (a dense table would have about 61
million entries at order 12).

## B-series

A B-series `a` is a vector over the standard numbering, `a[n] = a(tₙ)`.
The value `a(∅)` is not stored; it is 1 (central B-series). The series is

```math
B_h(a)\,y = a(\emptyset)\,y + \sum_{t} \frac{h^{|t|}\,a(t)}{\sigma(t)}\,F(t)(y),
```

the same normalisation as in BSeries.jl; the exact solution is
`E(t) = 1/t!`.

## Composition

[`compose`](@ref)`(a, b, T)` is the book's product `ab`: **first `a`, then
`b`** (Theorem 3.9),

```math
(ab)(t) = b(\emptyset)\,a(t) + \sum_{\emptyset < t' \le t} a(t \setminus t')\, b(t'),
```

where `t'` runs over the subtrees of `t` that contain the root, and
`a(F)` is the product of `a` over the trees of the forest `F`. This
agrees with `BSeries.compose(a, b)`.

The sum is taken from the table of [`splittings`](@ref), built once per
`Trees` object by the β-product recursion and reused for every
composition.

## Substitution

[`substitute`](@ref)`(b, a, T)` is `b ⋆ a` (argument order as in
`BSeries.substitute`): the method `a` applied to the vector field
`h⁻¹ B_h(b)`,

```math
(b \star a)(t) = \sum_{\text{partitions of } t} b(\text{piece}_1)\cdots b(\text{piece}_m)\; a(\text{skeleton}).
```

The partitions with their skeletons ([`tree_partitions`](@ref)) are built
by the same β-product recursion as the splittings.

## Number types

Functions are generic in the number type. The tableaux of the book are
given as rationals; [`exact_solution`](@ref) returns `Rational{BigInt}`,
and [`GeneralLinearMethod`](@ref) widens integer and rational
coefficients to `Rational{BigInt}`, so that no overflow occurs at high
orders. Floating-point tableaux stay floating point (use the `atol`
keyword of [`order_of_accuracy`](@ref)).
