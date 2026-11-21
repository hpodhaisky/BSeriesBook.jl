# Tutorial

This page is a guided tour of the package. Every code block is run when
the documentation is built, so the output shown is the output of the
current version.

```@example tut
using BSeriesBook
```

## Trees (Chapter 2)

`Trees(p)` builds all rooted trees up to order `p` in the book's standard
numbering (Algorithm 3 of the book), together with the tree functions
(Algorithm 4). Almost every other function takes this object as an
argument.

```@example tut
T = Trees(6)
```

A tree is its standard number `n`. Its B⁺ form, children, level sequence
and tree functions:

```@example tut
[(n, bracket(T, n), order(T, n), σ(T, n), γ(T, n)) for n in 1:9]
```

```@example tut
children(T, 6), tree_number(T, [1, 2]), level_sequence(T, 6)
```

The β-product `t_ℓ ∘ t_r` (graft `t_r` onto the root of `t_ℓ`) is a table
lookup:

```@example tut
n = beta_product(T, 2, 1)       # [τ] ∘ τ = [τ²]
bracket(T, n)
```

Counting trees (Algorithms 1 and 2):

```@example tut
count_trees(10)
```

```@example tut
U, UN = count_unrooted_trees(10)   # free trees, non-superfluous free trees
```

## B-series and composition (Chapter 3)

A B-series is a vector over the standard numbering, `a[n] = a(tₙ)`, with
`a(∅) = 1` implied. The exact solution is `E(t) = 1/t!`:

```@example tut
E = exact_solution(T)
E[1:8]
```

[`compose`](@ref)`(a, b, T)` is the book's product `ab`: first `a`, then
`b`. Two half steps of the exact flow give one full step:

```@example tut
half = [E[n] / 2^order(T, n) for n in 1:length(T)]
compose(half, half, T) == E
```

The inverse in the B-group:

```@example tut
a = [1 // (n + 1) for n in 1:length(T)]
compose(a, inverse(a, T), T) == identity_series(T)
```

Algorithm 6 gives the efficient representation λ̂ used in the book
to evaluate compositions:

```@example tut
λ̂ = lambda_hat(T, a, 4)
```

## Runge–Kutta methods (Chapters 1, 3, 5)

```@example tut
rk = kutta4()
rk.A, rk.b
```

Elementary weights Φ(t) and the order:

```@example tut
Φ = elementary_weights(rk, Trees(5))
Φ[1:8]'
```

```@example tut
[order_of_accuracy(m(), Trees(8)) for m in
    (euler, runge_midpoint, kutta4, kutta38, butcher6a)]
```

Any number type can be used, for example floating point:

```@example tut
g = gauss2()
order_of_accuracy(g, Trees(6); atol = 1e-12)
```

The stability function `R(z) = 1 + Σ bᵀ A^{k-1} 𝟏 zᵏ`:

```@example tut
stability_polynomial_coefficients(rk, 5)
```

## Modified equation and modifying integrator

The modified vector field `b` with `b ⋆ E = a`: the exact flow of
`y' = h⁻¹ B_h(b)(y)` reproduces the method. For the explicit Euler method,

```@example tut
modified_equation(euler(), Trees(4))'
```

i.e. `f − (h/2) f′f + h² ( f″(f,f)/12 + f′f′f/3 ) + …`. For the
classical fourth-order method the first nonzero correction appears at
order 5:

```@example tut
T7 = Trees(7)
b = modified_equation(kutta4(), T7)
n = findfirst(!iszero, b[2:end]) + 1
order(T7, n)
```

The modifying integrator (`b ⋆ a = E`) and the substitution law itself:

```@example tut
m = modifying_integrator(euler(), Trees(4))
substitute(m, elementary_weights(euler(), Trees(4)), Trees(4)) ==
    exact_solution(Trees(4))
```

## Conjugacy (Algorithm 7)

Two central B-series `a` and `b` are conjugate to order `p` if
`x a x⁻¹ = b + O(h^{p+1})` for some `x`. [`conjugate_order`](@ref) applies
Algorithm 7 ([`conjugacy_step`](@ref)) order by order. Conjugating the
exact solution with any `x` gives a series that is conjugate to `E` to
every order computed (here 6), although it differs from `E` already at
order 3:

```@example tut
T = Trees(6); E = exact_solution(T)
x = [order(T, n) == 2 ? 1//5 : 0//1 for n in 1:length(T)]
a = compose(compose(inverse(x, T), E, T), x, T)
first_difference = order(T, findfirst(n -> a[n] != E[n], 1:length(T)))
first_difference, conjugate_order(a, E, T)
```

The next page applies this to general linear methods.

## Working together with BSeries.jl

A B-series of BSeriesBook.jl can be converted to a `TruncatedBSeries` of
BSeries.jl via level sequences, for example to use functionality that only
BSeries.jl offers:

```julia
using BSeriesBook
import BSeries, RootedTrees

T = Trees(5)
a = elementary_weights(kutta4(), T)

s = BSeries.TruncatedBSeries{RootedTrees.RootedTree{Int, Vector{Int}},
                             eltype(a)}()
s[RootedTrees.rootedtree(Int[])] = one(eltype(a))
for n in 1:length(T)
    s[RootedTrees.rootedtree(level_sequence(T, n))] = a[n]
end

BSeries.order_of_accuracy(s)   # 4, as order_of_accuracy(kutta4(), T)
```
