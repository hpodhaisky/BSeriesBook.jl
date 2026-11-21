# [Underlying one-step method](@id glm)

This is the main addition of BSeriesBook.jl compared with BSeries.jl,
which has no support for multivalue methods.

## General linear methods

A general linear method `(A, U, B, V)` with `s` stages and `r` inputs
computes, in one step,

```math
Y = h\,A\,F(Y) + U\,y^{[n-1]}, \qquad y^{[n]} = h\,B\,F(Y) + V\,y^{[n-1]} .
```

Runge–Kutta methods (`r = 1`), linear multistep methods and many
multistage–multivalue methods are special cases (Book, Section 6.4).

## The underlying one-step method

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

## Example: an order-5 method of Chapter 6

Method (6.5 h), Case (ii), of the book (Section 6.5, "An example of
method derivation"), with four stages and two inputs:

```@example glm
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

```@example glm
T = Trees(7)
order_of_accuracy(glm, T)
```

The underlying one-step method agrees with the exact solution up to
order 5 and differs at order 6:

```@example glm
res = underlying_one_step(glm, T)
E = exact_solution(T)
[(p, all(res.E[n] == E[n] for n in trees_of_order(T, p))) for p in 1:7]
```

The starting method (in the original basis of the inputs) gives the
B-series of the second input quantity, here on the trees of order ≤ 4:

```@example glm
[res.ζ[n][2] for n in 1:8]'
```

## Checks

A Runge–Kutta method written as a general linear method must give its
own B-series:

```@example glm
rk = kutta4()
underlying_one_step(GeneralLinearMethod(rk), T).E == elementary_weights(rk, T)
```

A method that is not in the form (6.6 a), the third-order Adams–Bashforth
method written as a general linear method with four inputs (method
(6.4 m) of the book):

```@example glm
ab3 = GeneralLinearMethod(fill(0 // 1, 1, 1), [1 23//12 -4//3 5//12],
          reshape([0, 1, 0, 0] .// 1, 4, 1),
          [1 23//12 -4//3 5//12; 0 0 0 0; 0 1 0 0; 0 0 1 0])
is_partitioned_form(ab3), order_of_accuracy(ab3, Trees(6))
```

The stability matrix `M(z) = V + z B (I − z A)⁻¹ U`; its principal
eigenvalue is the stability function of `Ẽ`:

```@example glm
stability_matrix(glm, 1//10)
```

## Run time

For method (6.5 h)(ii) in `Rational{BigInt}`: [`underlying_one_step`](@ref)
takes about 5 ms at order 8 and 0.1 s at order 10;
`order_of_accuracy` about 0.4 s at order 8 and 4 s at order 10 (dominated
by [`conjugacy_step`](@ref), Algorithm 7).

## References

- J. C. Butcher, *B-Series*, Springer 2021, Sections 6.4–6.6.
- U. Kirchgraber, Multistep methods are essentially one-step methods,
  *Numer. Math.* 48 (1986), 85–90.
- D. Stoffer, General linear methods: connection to one step methods and
  invariant curves, *Numer. Math.* 64 (1993), 395–407.
