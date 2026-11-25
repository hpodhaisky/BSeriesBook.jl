# Chapters 1, 3 and 5: Runge–Kutta methods and their elementary weights.
#
# Book: Section 3.6 (elementary weights), Section 5.3 (stability function),
# Section 5.4 (explicit methods).

using LinearAlgebra: I, dot

"""
    RungeKuttaMethod(A, b, c = A * ones)

A Runge–Kutta method given by its tableau.  The entries may be rationals,
floating-point numbers or any other number type.
"""
struct RungeKuttaMethod{S}
    A::Matrix{S}
    b::Vector{S}
    c::Vector{S}
end

function RungeKuttaMethod(A::AbstractMatrix, b::AbstractVector,
                          c::AbstractVector = vec(sum(A; dims = 2)))
    S = promote_type(eltype(A), eltype(b), eltype(c))
    return RungeKuttaMethod{S}(Matrix{S}(A), Vector{S}(b), Vector{S}(c))
end

"""
    stages(rk), stages(glm)

Number of stages s of a Runge–Kutta or general linear method.
"""
stages(rk::RungeKuttaMethod) = length(rk.b)

"""
    stage_weights(rk, T) -> (ψ, φ)

`ψ[n]` and `φ[n] = A ψ[n]` for every tree t_n, from
    ψ(τ) = 𝟏,   ψ(t_ℓ ∘ t_r) = ψ(t_ℓ) ⊙ φ(t_r)   (componentwise).
`φ(t)` is the vector of internal (stage) weights of t.

Book: Section 3.6.
"""
function stage_weights(rk::RungeKuttaMethod, T::Trees)
    (; L, R) = T.table
    N = length(T)
    s = stages(rk)
    psi = Vector{Vector{eltype(rk.b)}}(undef, N)
    phi = similar(psi)
    psi[1] = ones(eltype(rk.b), s)
    phi[1] = rk.A * psi[1]
    for n in 2:N
        psi[n] = psi[L[n]] .* phi[R[n]]
        phi[n] = rk.A * psi[n]
    end
    return psi, phi
end

"""
    elementary_weights(rk, T) -> Vector

Φ(t_n) = bᵀ ψ(t_n) for every tree.  Book: Section 3.6.
"""
function elementary_weights(rk::RungeKuttaMethod, T::Trees)
    psi, _ = stage_weights(rk, T)
    return [dot(rk.b, psi[n]) for n in 1:length(T)]
end

"""
    bseries(rk, T) -> CentralBSeries

The B-series of one step of `rk`: a(∅) = 1, a(t) = Φ(t).
Book: Section 3.6.
"""
bseries(rk::RungeKuttaMethod, T::Trees) =
    CentralBSeries(T, elementary_weights(rk, T))

"""
    order_of_accuracy(rk::RungeKuttaMethod, T; atol = 0) -> Int

The order of `rk`, at most the order of `T`: the order of accuracy of
its B-series.  Use `atol > 0` for floating-point tableaux.

Book: Section 3.6 (order conditions).
"""
order_of_accuracy(rk::RungeKuttaMethod, T::Trees; atol = 0) =
    order_of_accuracy(bseries(rk, T); atol)

"""
    stability_polynomial_coefficients(rk, m) -> Vector

The coefficients r₀, r₁, …, r_m of R(z) = 1 + z bᵀ(I - zA)⁻¹𝟏
= Σ_k bᵀA^{k-1}𝟏 z^k (r₀ = 1).  For an explicit s-stage method R is a
polynomial of degree ≤ s, and m = s gives all of it.

Book: Section 5.3.
"""
function stability_polynomial_coefficients(rk::RungeKuttaMethod, m::Integer)
    v = ones(eltype(rk.b), stages(rk))
    r = [one(eltype(rk.b))]
    for _ in 1:m
        push!(r, dot(rk.b, v))
        v = rk.A * v
    end
    return r
end

# ---------------------------------------------------------------------------
# Tableaux used in the book
# ---------------------------------------------------------------------------

"""Explicit Euler method."""
euler() = RungeKuttaMethod(fill(0 // 1, 1, 1), [1 // 1])

"""Runge's second-order (midpoint) method.  Book: Section 1.5."""
runge_midpoint() = RungeKuttaMethod([0//1 0; 1//2 0], [0 // 1, 1])

"""The classical fourth-order method of Kutta.  Book: Sections 1.5, 5.4."""
kutta4() = RungeKuttaMethod(
    [0//1 0 0 0
     1//2 0 0 0
     0 1//2 0 0
     0 0 1 0],
    [1 // 6, 1 // 3, 1 // 3, 1 // 6])

"""Kutta's 3/8-rule, fourth order.  Book: Section 5.4."""
kutta38() = RungeKuttaMethod(
    [0//1 0 0 0
     1//3 0 0 0
     -1//3 1 0 0
     1 -1 1 0],
    [1 // 8, 3 // 8, 3 // 8, 1 // 8])

"""
    gauss2(S = Float64)

The two-stage Gauss method of order 4 (coefficients involve √3).
Book: Section 5.6.
"""
function gauss2(::Type{S} = Float64) where {S}
    r = sqrt(S(3)) / 6
    A = [S(1) / 4  S(1) / 4 - r
         S(1) / 4 + r  S(1) / 4]
    return RungeKuttaMethod(A, [S(1) / 2, S(1) / 2])
end

"""
    butcher6a()

A seven-stage method of order 6 with c₄ = c₃/(2 − 10c₃ + 15c₃²),
proposed for the 2nd edition (J. C. Butcher, "Efficient order 6
methods").
"""
butcher6a() = RungeKuttaMethod(
    [0//1 0 0 0 0 0 0
     1//3 0 0 0 0 0 0
     0 2//3 0 0 0 0 0
     1//12 1//3 -1//12 0 0 0 0
     -1//16 9//8 -3//16 -3//8 0 0 0
     0 9//8 -3//8 -3//4 1//2 0 0
     9//44 -9//11 63//44 18//11 0 -16//11 0],
    [11 // 120, 0, 27 // 40, 27 // 40, -4 // 15, -4 // 15, 11 // 120])
