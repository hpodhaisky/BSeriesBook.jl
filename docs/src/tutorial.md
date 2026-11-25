# Tutorial

This page explains the package step by step, following the chapters of
the book: first the conventions — how trees and B-series are
represented — and then the main computations. Every code block is run
when the documentation is built, so the output shown is the output of
the current version. The [Reference](@ref) lists all functions.

```@example tut
using BSeriesBook
```

## [Trees and the standard numbering (Chapter 2)](@id standard-numbering)

The book numbers the rooted trees by Algorithm 3 ("Generate trees",
Section 2.5). Every tree except `t₁ = τ` is the β-product
`tₙ = t_{L[n]} ∘ t_{R[n]}`: the tree `t_{R[n]}` is grafted onto the root
of `t_{L[n]}`, and `R[n] ≤ R[L[n]]` makes the decomposition unique. Within
each order the trees are numbered in the order in which the algorithm
produces them.

`Trees(p)` builds all trees up to order `p` in this numbering, together
with the tree functions (Algorithm 4). Almost every other function takes
this object as an argument.

```@example tut
T = Trees(6)
```

A tree is a [`Tree`](@ref): `T[n]` wraps the standard number `n`
together with its `Trees`. It is shown with its number and its B⁺ form:

```@example tut
t = T[17]
```

Functions on trees take a `Tree`; [`number`](@ref) gives back `n`. The
type costs nothing at run time: a loop over `Tree`s compiles to the same
code as a loop over the numbers.

```@example tut
[(t, order(t), symmetry(t), factorial(t)) for t in T if order(t) <= 4]
```

```@setup trees
using BSeriesBook

# Positions of the vertices of tree n: leaves get consecutive x positions,
# an inner vertex sits above the middle of its children.
function place!(pos, edges, T, n, x0, depth)
    kids = [number(k) for k in children(T[n])]
    if isempty(kids)
        push!(pos, (x0 + 0.5, depth))
        return length(pos), x0 + 1
    end
    push!(pos, (0.0, depth)); me = length(pos)
    xs = Float64[]; x = x0
    for k in kids
        c, x = place!(pos, edges, T, k, x, depth + 1)
        push!(edges, (me, c)); push!(xs, pos[c][1])
    end
    pos[me] = ((first(xs) + last(xs)) / 2, depth)
    return me, x
end

function tree_svg(T, n; dx = 13, dy = 16, r = 3.2, pad = 5)
    pos = Tuple{Float64,Int}[]; edges = Tuple{Int,Int}[]
    _, w = place!(pos, edges, T, n, 0.0, 0)
    h = maximum(last, pos)
    W = w * dx + 2pad - dx + 2r; H = h * dy + 2pad + 2r
    X(p) = pad + r + (p[1] - 0.5) * dx
    Y(p) = H - pad - r - p[2] * dy
    s = IOBuffer()
    print(s, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$(W)\" ",
          "height=\"$(H)\" style=\"vertical-align:middle\">")
    for (a, b) in edges
        print(s, "<line x1=\"$(X(pos[a]))\" y1=\"$(Y(pos[a]))\" ",
              "x2=\"$(X(pos[b]))\" y2=\"$(Y(pos[b]))\" ",
              "stroke=\"currentColor\" stroke-width=\"1.3\"/>")
    end
    for p in pos
        print(s, "<circle cx=\"$(X(p))\" cy=\"$(Y(p))\" r=\"$r\" ",
              "fill=\"currentColor\"/>")
    end
    print(s, "</svg>")
    return String(take!(s))
end

function html_table(header, rows)
    s = IOBuffer()
    print(s, "<table><thead><tr>")
    foreach(h -> print(s, "<th style=\"text-align:center\">", h, "</th>"),
            header)
    print(s, "</tr></thead><tbody>")
    for row in rows
        print(s, "<tr>")
        foreach(c -> print(s, "<td style=\"text-align:center\">", c,
                           "</td>"), row)
        print(s, "</tr>")
    end
    print(s, "</tbody></table>")
    return HTML(String(take!(s)))
end

T = Trees(5)
L, R = T.table.L, T.table.R
lr(n, v) = n == 1 ? "–" : "t<sub>$(v[n])</sub>"
```


The trees up to order 5 with their tree functions (the drawings have the
root at the bottom, as in the book):

| Order | 1 | 2 | 3 | 4 | 5 |
|:------|:-:|:-:|:-:|:-:|:-:|
| number of trees | 1 | 1 | 2 | 4 | 9 |
| standard numbers | 1 | 2 | 3–4 | 5–8 | 9–17 |

```@example trees
html_table(["n", "tree", "B⁺ form", "|t|", "L[n]", "R[n]",
            "σ(t)", "t!", "α(t)", "β(t)"],
    [["t<sub>$n</sub>", tree_svg(T, n), bracket(T[n]), order(T[n]),
      lr(n, L), lr(n, R), symmetry(T[n]), factorial(T[n]), alpha(T[n]), beta(T[n])]
     for n in 1:length(T)])
```

Here σ(t) is the symmetry, t! the tree factorial,
α(t) = |t|!/(σ(t) t!) the number of essentially different labellings
and β(t) = |t|!/σ(t) (Book, Section 2.5).

### Product table

The β-product `t_ℓ ∘ t_r` for all pairs with |t_ℓ| + |t_r| ≤ 5 (row `ℓ`,
column `r`), the table `prod` of Algorithm 3:

```@example trees
html_table(vcat(["l \\ r"], ["t<sub>$r</sub>" for r in 1:8]),
    [vcat(["t<sub>$l</sub>"],
          [order(T[l]) + order(T[r]) <= 5 ?
               "t<sub>$(number(beta_product(T[l], T[r])))</sub>" : ""
           for r in 1:8])
     for l in 1:8])
```

The β-product `s ∘ t` is a lookup in this table:

```@example tut
beta_product(T[2], T[1])       # [τ] ∘ τ = [τ²]
```

A forest is a [`Forest`](@ref), the sorted vector of its trees, shown as
a monomial. The children B⁻(t) of `t = [t₁ ⋯ t_m]` form a forest, and
[`bplus`](@ref) is the inverse:

```@example tut
children(T[6]), bplus(T, Forest([T[1], T[2]]))
```

The splittings (t′, t ∖ t′) of a tree — a subtree t′ containing the root
and the forest of the removed parts — are needed for the composition:

```@example tut
splittings(T[4])
```

Counting trees (Algorithms 1 and 2):

```@example tut
count_trees(10)
```

```@example tut
U, UN = BSeriesBook.count_unrooted_trees(10)   # free trees, non-superfluous free trees
```

## B-series (Chapter 3)

A B-series is a mapping a: T# → ℝ, T# = {∅} ∪ T, and stands for the
series

```math
B_h(a)\,y = a(\emptyset)\,y + \sum_{t} \frac{h^{|t|}\,a(t)}{\sigma(t)}\,F(t)(y),
```

the same normalisation as in BSeries.jl. In the package it is a vector
over the standard numbering, `a[t] = a(t)`, truncated at the order of its
`Trees`; the value on the empty tree is `a[empty_tree]` (or `a[∅]`). The
type records a(∅), following the sets 𝖡\*, 𝖡 and 𝖡⁰ of Section 3.4 of the
book: [`GeneralBSeries`](@ref), [`CentralBSeries`](@ref) (a(∅) = 1) and
[`VectorFieldBSeries`](@ref) (a(∅) = 0).

The exact solution `E(t) = 1/t!` is a central B-series:

```@example tut
E = exact_solution(T)
```

```@example tut
E[empty_tree], E[T[3]], coefficients(E)[1:4]
```

The type is decided by the structure of a computation, never by the
value of a(∅): the difference of two central series lies in 𝖡⁰, the sum
of two central series in 𝖡\*. A conversion such as `CentralBSeries(g)`
checks the value of a(∅). Equality `==` and `≈` compare a(∅) as well;
indexing with a range, broadcasting and other generic vector operations
act on the coefficients a(tₙ) only and return plain vectors.

```@example tut
a = bseries(kutta4(), T)
typeof(a - E), typeof(a + E), (a + E)[empty_tree]
```

## Composition

[`compose`](@ref)`(a, b)`, also written `a * b`, is the book's product
`ab`: **first `a`, then `b`** (Theorem 3.9), a map 𝖡 × 𝖡\* → 𝖡\*,

```math
(ab)(\emptyset) = b(\emptyset), \qquad
(ab)(t) = b(\emptyset)\,a(t) + \sum_{\emptyset < t' \le t} a(t \setminus t')\, b(t'),
```

where `t'` runs over the subtrees of `t` that contain the root (the
splittings above), and `a(F)` is the product of `a` over the trees of
the forest `F`. This agrees with `BSeries.compose(a, b)`. The first
factor must be central (otherwise `B_h(b) ∘ B_h(a)` is not a B-series);
the result has the type of `b`. Two half steps of the exact flow give
one full step:

```@example tut
half = rescale(E, 1//2)          # E^(1/2)(t) = 2^-|t| / t!
half * half == E
```

The B-group: inverse, identity and powers.

```@example tut
a * inv(a) == one(a), a^3 == a * a * a, a^-1 == inverse(a)
```

A product with a non-central second factor has the type of the second
factor:

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

The B-series of `hf(y₀)` is D ∈ 𝖡⁰ (D(τ) = 1, D(t) = 0 otherwise), and
for a central `y` the series of `hf(B_h(y)y₀)` is `yD`, with
`(yD)([t₁ ⋯ t_m]) = y(t₁) ⋯ y(t_m)`. [`slope`](@ref)`(y)` computes `yD`;
for the identity it is D itself:

```@example tut
D = slope(identity_series(T))
coefficients(D)[1:4]
```

`slope(y, gamma)` solves the implicit relation `w = y + γ hf(w)` in one
pass over the standard numbering: it returns `dw = wD` with
`w = y + γ dw`. With `y` the identity and `γ = 1`, `w` is the implicit
Euler method:

```@example tut
w = identity_series(T) + slope(identity_series(T), 1)
w == bseries(RungeKuttaMethod(fill(1//1, 1, 1), [1//1]), T)
```

## Runge–Kutta methods (Chapters 1, 3, 5)

```@example tut
rk = kutta4()
rk.A, rk.b
```

Elementary weights Φ(t), the B-series of one step and the order:

```@example tut
Phi = elementary_weights(rk, Trees(5))
Phi[1:8]'
```

```@example tut
[order_of_accuracy(m(), Trees(8)) for m in
    (euler, runge_midpoint, kutta4, kutta38, butcher6a)]
```

A diagonally implicit method, stage by stage with `slope`:
`Fᵢ = slope(yᵢ, aᵢᵢ)` with `yᵢ = 1 + Σ_{j<i} aᵢⱼ Fⱼ`:

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

The stability function `R(z) = 1 + Σ bᵀ A^{k-1} 𝟏 zᵏ`:

```@example tut
stability_polynomial_coefficients(rk, 5)
```

## Number types

Functions are generic in the number type. The tableaux of the book are
given as rationals; [`exact_solution`](@ref) returns `Rational{BigInt}`,
and [`GeneralLinearMethod`](@ref) widens integer and rational
coefficients to `Rational{BigInt}`, so that no overflow occurs at high
orders. Floating-point tableaux stay floating point; use the `atol`
keyword of [`order_of_accuracy`](@ref):

```@example tut
order_of_accuracy(gauss2(), Trees(6); atol = 1e-12)
```

## Substitution and the modified equation

[`substitute`](@ref)`(b, a)` is `b ⋆ a` for `b` in 𝖡⁰ (argument order as
in `BSeries.substitute`): the method `a` applied to the vector field
`h⁻¹ B_h(b)`,

```math
(b \star a)(t) = \sum_{\text{partitions of } t} b(\text{piece}_1)\cdots b(\text{piece}_m)\; a(\text{skeleton}),
```

with `(b ⋆ a)(∅) = a(∅)`; the result has the type of `a`. A partition
cuts a set of edges of `t`; the skeleton is `t` with every piece
contracted to a vertex ([`tree_partitions`](@ref)):

```@example tut
tree_partitions(T[4])
```

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
bm = modified_equation(kutta4(), T7)
order(first(t for t in T7 if order(t) > 1 && bm[t] != 0))
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
x = CentralBSeries(t -> order(t) == 2 ? 1//5 : 0//1, T)
ax = inv(x) * E * x
first_difference = order(first(t for t in T if ax[t] != E[t]))
first_difference, conjugate_order(ax, E)
```

## [General linear methods and the underlying one-step method (Chapter 6)](@id glm)

This is the main addition of BSeriesBook.jl compared with BSeries.jl,
which has no support for multivalue methods.

A general linear method `(A, U, B, V)` with `s` stages and `r` inputs
computes, in one step,

```math
Y = h\,A\,F(Y) + U\,y^{[n-1]}, \qquad y^{[n]} = h\,B\,F(Y) + V\,y^{[n-1]} .
```

Runge–Kutta methods (`r = 1`), linear multistep methods and many
multistage–multivalue methods are special cases (Book, Section 6.4).

### The underlying one-step method

The `r` input quantities approximate `r` different B-series of the
solution. Kirchgraber (1986) and Stoffer (1993) showed that there is a
**starting method** `S` and a **one-step method** `Ẽ` (both formal
B-series) with

```math
M \circ S = S \circ \tilde E ,
```

where `M` is one step of the general linear method: starting from `S(y₀)`,
the method behaves exactly like the one-step method `Ẽ`, observed through
`S`. The **order** of the general linear method is the order to which `Ẽ`
is *conjugate* to the exact solution `E` (Book, Section 6.6).

For a method in the partitioned form (6.6 a) of the book,

```math
\begin{bmatrix} A & U \\ B & V \end{bmatrix} =
\begin{bmatrix} A & \mathbf 1 & \hat U \\ b^T & 1 & 0^T \\ \hat B & 0 & \hat V \end{bmatrix},
\qquad 1 \notin \sigma(\hat V),
```

`Ẽ` and the starting method `(identity, ξ)` are computed tree by tree
from the book's recursion (6.6 d) (details in Table 18 of the book):

```math
\begin{aligned}
\eta(t) &= A\,(\eta D)(t) + \hat U\,\xi(t), \\
\tilde E(t) &= b^T (\eta D)(t), \\
\xi(t) &= (I - \hat V)^{-1}\Bigl(\hat B\,(\eta D)(t) - \sum_{t' < t} \tilde E(t \setminus t')\,\xi(t')\Bigr),
\end{aligned}
```

with `(ηD)(t_ℓ ∘ t_r) = (ηD)(t_ℓ) ⊙ η(t_r)`. The sum is taken from the
same table of splittings as [`compose`](@ref). Every preconsistent method
whose matrix `V` has 1 as a simple eigenvalue is first brought into the
form (6.6 a) by a change of basis of the input vector
([`partitioned_form`](@ref)).

### Example: an order-5 method of Chapter 6

Method (6.5 h), Case (ii), of the book (Section 6.5, "An example of
method derivation"), with four stages and two inputs:

```@example tut
using BSeriesBook

A = [0 0 0 0; 3//4 0 0 0; 93//250 -9//125 0 0;
     8881//8640 1069//3240 -1855//5184 0]
U = [1 0; 1 1; 1 44//125; 1 691//540]
B = [5//54 32//81 250//567 1//14;
     -19//16 -37//6 1175//336 27//7]
V = [1 0; 0 0]
glm = GeneralLinearMethod(A, U, B, V)
stages(glm), inputs(glm), is_partitioned_form(glm)
```

```@example tut
T = Trees(7)
order_of_accuracy(glm, T)
```

The underlying one-step method agrees with the exact solution up to
order 5 and differs at order 6:

```@example tut
res = underlying_one_step(glm, T)
E = exact_solution(T)
[(p, all(res.E[t] == E[t] for t in trees_of_order(T, p))) for p in 1:7]
```

The starting method (in the original basis of the inputs) consists of
one B-series per input, `res.S[i]`, with `res.S[i][empty_tree] = u[i]`; here the
second input quantity:

```@example tut
res.S[2]
```

`M ∘ S = S ∘ Ẽ` can be read off component by component: one step of the
underlying method, followed by the starting method, is
`res.E * res.S[i]`.

### Checks

A Runge–Kutta method written as a general linear method must give its
own B-series:

```@example tut
rk = kutta4()
underlying_one_step(GeneralLinearMethod(rk), T).E == bseries(rk, T)
```

A method that is not in the form (6.6 a), the third-order Adams–Bashforth
method written as a general linear method with four inputs (method
(6.4 m) of the book):

```@example tut
ab3 = GeneralLinearMethod(fill(0 // 1, 1, 1), [1 23//12 -4//3 5//12],
          reshape([0, 1, 0, 0] .// 1, 4, 1),
          [1 23//12 -4//3 5//12; 0 0 0 0; 0 1 0 0; 0 0 1 0])
is_partitioned_form(ab3), order_of_accuracy(ab3, Trees(6))
```

The stability matrix `M(z) = V + z B (I − z A)⁻¹ U`; its principal
eigenvalue is the stability function of `Ẽ`:

```@example tut
BSeriesBook.stability_matrix(glm, 1//10)
```

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
    s[RootedTrees.rootedtree(BSeriesBook.level_sequence(t))] = a[t]
end

BSeries.order_of_accuracy(s)   # 4, as order_of_accuracy(kutta4(), T)
```
