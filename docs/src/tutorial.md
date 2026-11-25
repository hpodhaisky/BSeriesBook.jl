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

A tree is a [`Tree`](@ref): `T[n]` is the tree with standard number
`n`. It is shown with its number and its B⁺ form:

```@example tut
t = T[17]
```

Tree functions take the tree, not a number:

```@example tut
[(t, order(t), symmetry(t), factorial(t)) for t in T if order(t) <= 4]
```

Children, B⁺ of a forest, level sequence and the standard number:

```@example tut
children(T[6]), bplus(T, [T[1], T[2]]), level_sequence(T[6]), number(T[6])
```

The β-product `s ∘ t` (graft `t` onto the root of `s`) is a table
lookup:

```@example tut
beta_product(T[2], T[1])       # [τ] ∘ τ = [τ²]
```

Counting trees (Algorithms 1 and 2):

```@example tut
count_trees(10)
```

```@example tut
U, UN = count_unrooted_trees(10)   # free trees, non-superfluous free trees
```

## B-series (Chapter 3)

A B-series is a vector over the standard numbering, `a[n] = a(tₙ)`,
whose type records the value on the empty tree, following the sets 𝖡\*,
𝖡 and 𝖡⁰ of the book (see [Conventions](@ref)). The exact solution
`E(t) = 1/t!` is a central B-series:

```@example tut
E = exact_solution(T)
```

```@example tut
E[empty_tree], E[3], coefficients(E)[1:4]
```

The type of a linear combination follows from a(∅): the difference of
two central series lies in 𝖡⁰, a sum of two central series in 𝖡\*.

```@example tut
a = bseries(kutta4(), T)
typeof(a - E), typeof(a + E), (a + E)[empty_tree]
```

## Composition

[`compose`](@ref)`(a, b)`, or `a * b`, is the book's product `ab`: first
`a`, then `b`. The first factor must be central. Two half steps of the
exact flow give one full step:

```@example tut
half = rescale(E, 1//2)          # E^(1/2)(t) = 2^-|t| / t!
half * half == E
```

The B-group: inverse, identity and powers.

```@example tut
a * inv(a) == one(a), a^3 == a * a * a, a^-1 == inverse(a)
```

The product is also defined when the second factor is not central; it
has the type of the second factor:

```@example tut
g = GeneralBSeries(t -> 1//number(t), T, 2)
typeof(a * g), (a * g)[empty_tree]
```

Algorithm 6 gives the efficient representation λ̂ used in the book
to evaluate compositions:

```@example tut
lhat = lambda_hat(a, 4)
```

## The map slope and implicit stages

`slope(y)` is the series `yD` of `hf(B_h(y)y₀)`; for the identity it is
D itself:

```@example tut
D = slope(identity_series(T))
coefficients(D)[1:4]
```

`slope(y, aii)` solves the implicit relation of a diagonally implicit
stage, `dw = wD` with `w = y + aii·dw`. With `y` the identity and
`aii = 1` the stage `w` is the implicit Euler method:

```@example tut
w = identity_series(T) + slope(identity_series(T), 1)
w == bseries(RungeKuttaMethod(fill(1//1, 1, 1), [1//1]), T)
```

A diagonally implicit method stage by stage, `Fᵢ = slope(yᵢ, aᵢᵢ)` with
`yᵢ = 1 + Σ_{j<i} aᵢⱼ Fⱼ`:

```@example tut
A = [1//3 0 0; 1//2 1//4 0; -1//5 2//3 1//2]; b = [1//4, 1//2, 1//4]
id = identity_series(T)
F = VectorFieldBSeries[]
for i in 1:3
    y = id + sum((A[i, j] * F[j] for j in 1:i-1); init = zero(id))
    push!(F, slope(y, A[i, i]))
end
id + sum(b[i] * F[i] for i in 1:3) == bseries(RungeKuttaMethod(A, b), T)
```

## Runge–Kutta methods (Chapters 1, 3, 5)

```@example tut
rk = kutta4()
rk.A, rk.b
```

Elementary weights Φ(t), the B-series of one step and the order:

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

The modified vector field `b` ∈ 𝖡⁰ with `b ⋆ E = a`: the exact flow of
`y' = h⁻¹ B_h(b)(y)` reproduces the method. For the explicit Euler method,

```@example tut
modified_equation(euler(), Trees(4))
```

i.e. `f − (h/2) f′f + h² ( f″(f,f)/12 + f′f′f/3 ) + …`. For the
classical fourth-order method the first nonzero correction appears at
order 5:

```@example tut
T7 = Trees(7)
b = modified_equation(kutta4(), T7)
t = first(t for t in T7 if order(t) > 1 && b[t] != 0)
order(t)
```

The modifying integrator (`b ⋆ a = E`) and the substitution law itself:

```@example tut
T4 = Trees(4)
m = modifying_integrator(euler(), T4)
substitute(m, bseries(euler(), T4)) == exact_solution(T4)
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
x = CentralBSeries(t -> order(t) == 2 ? 1//5 : 0//1, T)
a = inv(x) * E * x
first_difference = order(first(t for t in T if a[t] != E[t]))
first_difference, conjugate_order(a, E)
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
a = bseries(kutta4(), T)

s = BSeries.TruncatedBSeries{RootedTrees.RootedTree{Int, Vector{Int}},
                             eltype(a)}()
s[RootedTrees.rootedtree(Int[])] = a[empty_tree]
for t in T
    s[RootedTrees.rootedtree(level_sequence(t))] = a[t]
end

BSeries.order_of_accuracy(s)   # 4, as order_of_accuracy(kutta4(), T)
```
