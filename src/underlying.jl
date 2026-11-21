# Chapter 6: general linear methods and the underlying one-step method.
#
# A general linear method (A, U, B, V) with s stages and r inputs computes
#     Y = h A F(Y) + U y^[n-1],      y^[n] = h B F(Y) + V y^[n-1].
# For a method in the partitioned form (6.6 a) of the book,
#     [A U; B V] = [A  𝟏  Û;  bᵀ  1  0ᵀ;  B̂  0  V̂],    1 ∉ σ(V̂),
# the underlying one-step method Ẽ and the starting method (identity, ξ)
# satisfy M ∘ S = S ∘ Ẽ as formal series; they are computed tree by tree
# from (6.6 d):
#     η(t) = A (ηD)(t) + Û ξ(t),
#     Ẽ(t) = bᵀ (ηD)(t),
#     ξ(t) = (I - V̂)⁻¹ ( B̂ (ηD)(t) - Σ_{t'<t} Ẽ(t ∖ t') ξ(t') ),
# with (ηD)([t₁ ⋯ t_m]) = η(t₁) ⋯ η(t_m), η(∅) = 𝟏, ξ(∅) = 0, Ẽ(∅) = 1.
# Every preconsistent method whose matrix V has 1 as a simple eigenvalue
# is brought into this form by a change of basis of the input vector.
#
# Book: Section 6.6, (6.6 a), (6.6 d) and Table 18; Kirchgraber (1986),
# Stoffer (1993).

using LinearAlgebra: lu

"""
    GeneralLinearMethod(A, U, B, V)

A general linear method with `s = size(A, 1)` stages and `r = size(V, 1)`
inputs.  Book: Section 6.4.
"""
struct GeneralLinearMethod{S}
    A::Matrix{S}
    U::Matrix{S}
    B::Matrix{S}
    V::Matrix{S}
end

function GeneralLinearMethod(A::AbstractMatrix, U::AbstractMatrix,
                             B::AbstractMatrix, V::AbstractMatrix)
    S = exact(promote_type(eltype(A), eltype(U), eltype(B), eltype(V)))
    return GeneralLinearMethod{S}(Matrix{S}(A), Matrix{S}(U),
                                  Matrix{S}(B), Matrix{S}(V))
end

# Rationals are widened to BigInt to avoid overflow in the B-series.
exact(::Type{<:Union{Integer,Rational}}) = Rational{BigInt}
exact(S::Type) = S

"""
    GeneralLinearMethod(rk::RungeKuttaMethod)

A Runge–Kutta method as a general linear method with r = 1.
"""
GeneralLinearMethod(rk::RungeKuttaMethod) =
    GeneralLinearMethod(rk.A, ones(eltype(rk.b), stages(rk), 1),
                        reshape(rk.b, 1, :), ones(eltype(rk.b), 1, 1))

stages(glm::GeneralLinearMethod) = size(glm.A, 1)

"""
    inputs(glm)

Number r of input (and output) quantities of a general linear method.
"""
inputs(glm::GeneralLinearMethod) = size(glm.V, 1)

"""
    stability_matrix(glm, z) -> Matrix

M(z) = V + z B (I - z A)⁻¹ U, the matrix of one step for y′ = λy, z = hλ.
Book: Section 6.4.
"""
function stability_matrix(glm::GeneralLinearMethod, z)
    (; A, U, B, V) = glm
    s = stages(glm)
    Id = [i == j ? one(z) : zero(z) for i in 1:s, j in 1:s]
    return V + z * B * ((Id - z * A) \ U)
end

"""
    is_partitioned_form(glm) -> Bool

Whether `glm` has the partitioned form (6.6 a) of the book: U = [𝟏 Û], B = [bᵀ; B̂],
V = [1 0ᵀ; 0 V̂], with 1 not an eigenvalue of V̂.
"""
function is_partitioned_form(glm::GeneralLinearMethod)
    (; U, B, V) = glm
    r = inputs(glm)
    all(==(1), U[:, 1]) || return false
    V[1, 1] == 1 || return false
    all(iszero, V[1, 2:r]) && all(iszero, V[2:r, 1]) || return false
    r == 1 && return true
    Vh = V[2:r, 2:r]
    Id = [i == j ? one(eltype(V)) : zero(eltype(V)) for i in 1:r-1, j in 1:r-1]
    return pivot_columns(Id - Vh) == 1:r-1
end

"""
    partitioned_form(glm) -> (glm′, Q)

An equivalent method in the partitioned form (6.6 a), with inputs ỹ = Q⁻¹ y:
A′ = A, U′ = U Q, B′ = Q⁻¹ B, V′ = Q⁻¹ V Q.  The first column of Q is
the preconsistency vector u (V u = u, U u = 𝟏); the others span the
range of V - I, the invariant subspace for the other eigenvalues of V.
The first component of ỹ is wᵀ y, w the left eigenvector of V for the
eigenvalue 1 with wᵀu = 1.  Requires 1 to be a simple eigenvalue of V.

Book: Section 6.6, (6.6 a).
"""
function partitioned_form(glm::GeneralLinearMethod)
    (; A, U, B, V) = glm
    r, s = inputs(glm), stages(glm)
    S = eltype(V)
    Id = [i == j ? one(S) : zero(S) for i in 1:r, j in 1:r]
    # preconsistency vector: (V - I) u = 0, U u = 𝟏
    tol = S <: AbstractFloat ? sqrt(eps(S)) : 0
    u = solve_exact([V - Id; U], [zeros(S, r); ones(S, s)]; atol = tol)
    u === nothing && error("method is not preconsistent")
    cols = pivot_columns(V - Id)
    length(cols) == r - 1 ||
        error("1 is not a simple eigenvalue of V (geometric multiplicity)")
    Q = [u (V - Id)[:, cols]]
    length(pivot_columns(Q)) == r ||
        error("1 is not a simple eigenvalue of V (algebraic multiplicity)")
    glm′ = GeneralLinearMethod(A, U * Q, Q \ B, Q \ (V * Q))
    return glm′, Q
end

# Indices of the pivot columns of M (Gaussian elimination, exact for
# rationals; for floating point a column is a pivot if its entry is
# larger than a small tolerance).
function pivot_columns(M::AbstractMatrix)
    A = copy(M)
    m, n = size(A)
    tol = eltype(A) <: AbstractFloat ? 1e3 * eps(eltype(A)) * max(1, maximum(abs, A; init = 0)) : 0
    pivots = Int[]
    row = 1
    for col in 1:n
        row > m && break
        _, k = findmax(i -> abs(A[i, col]), row:m)
        k += row - 1                                  # row of the pivot
        abs(A[k, col]) > tol || continue
        A[[row, k], :] = A[[k, row], :]
        for i in row+1:m
            A[i, :] .-= (A[i, col] / A[row, col]) .* A[row, :]
        end
        push!(pivots, col)
        row += 1
    end
    return pivots
end

"""
    underlying_one_step(glm, T) -> (; E, ξ, η, ηD, ζ)

The underlying one-step method of `glm` and the starting method that
goes with it, to the order of `T`:
- `E[n]` = Ẽ(t_n), the B-series of the underlying one-step method
  (Ẽ(∅) = 1, not stored);
- `ξ[n]` = ξ(t_n) ∈ ℝ^{r-1}, the non-principal part of the starting
  method of the partitioned form (ξ(∅) = 0);
- `η[n]`, `ηD[n]` ∈ ℝ^s, the stages and stage derivatives;
- `ζ[n]` ∈ ℝ^r, the starting method in the original input basis,
  ζ(∅) = u, ζ(t) = Q [0; ξ(t)]: its principal part wᵀζ is the identity.
A method that is not in the partitioned form (6.6 a) is first
transformed with [`partitioned_form`](@ref).

Book: Section 6.6, recursion (6.6 d) and Table 18.
"""
function underlying_one_step(glm::GeneralLinearMethod, T::Trees)
    glm′, Q = is_partitioned_form(glm) ?
        (glm, [i == j ? one(eltype(glm.V)) : zero(eltype(glm.V))
               for i in 1:inputs(glm), j in 1:inputs(glm)]) :
        partitioned_form(glm)
    E, ξ, η, ηD = underlying_recursion(glm′, T)
    u = Q[:, 1]
    ζ = [Q * [zero(eltype(u)); ξ[n]] for n in 1:length(T)]
    return (; E, ξ, η, ηD, ζ, u, Q)
end

# Recursion (6.6 d) for a method in the partitioned form (6.6 a).
function underlying_recursion(glm::GeneralLinearMethod, T::Trees)
    (; A, U, B, V) = glm
    (; L, R) = T.table
    r, s = inputs(glm), stages(glm)
    S = eltype(V)
    b = B[1, :]
    Uh, Bh, Vh = U[:, 2:r], B[2:r, :], V[2:r, 2:r]
    Id = [i == j ? one(S) : zero(S) for i in 1:r-1, j in 1:r-1]
    IminusV = r > 1 ? lu(Id - Vh) : nothing
    solve(v) = r > 1 ? IminusV \ v : v
    N = length(T)
    Sp = splittings(T)
    E = zeros(S, N)
    ξ = Vector{Vector{S}}(undef, N)
    η = Vector{Vector{S}}(undef, N)
    ηD = Vector{Vector{S}}(undef, N)
    for n in 1:N
        # [book:MV28:begin]
        ηD[n] = n == 1 ? ones(S, s) : ηD[L[n]] .* η[R[n]]
        E[n] = b' * ηD[n]
        w = Bh * ηD[n]
        for (k, F, m) in Sp[n]
            k == n && continue
            w -= m * forest_value(E, F) * ξ[k]
        end
        ξ[n] = solve(w)
        η[n] = A * ηD[n] + Uh * ξ[n]
        # [book:MV28:end]
    end
    return E, ξ, η, ηD
end

"""
    conjugate_order(a, b, T) -> p

The largest p ≤ p_max such that a is conjugate to b to order p, i.e.
x a x⁻¹ = b + O_{p+1} for some x, found by repeated use of Algorithm 7
of the book, [`conjugacy_step`](@ref) (one order at a time, free
parameters set to 0).  For a = Ẽ, the
underlying one-step method of a general linear method, and b = E, this
is the order of the method (with a suitable starting method).

The test is exact: use rational (or other exact) coefficients.

Book: Section 6.6 ("Test for conjugacy").
"""
function conjugate_order(a::AbstractVector, b::AbstractVector, T::Trees)
    a[1] == b[1] || return 0
    a′ = a
    for p in 2:T.table.p_max
        res = conjugacy_step(a′, b, T, p)
        res[1] || return p - 1
        a′ = res[2]
    end
    return T.table.p_max
end

"""
    order_of_accuracy(glm::GeneralLinearMethod, T)

The order of a general linear method: the conjugate order of its
underlying one-step method to the exact solution.
"""
order_of_accuracy(glm::GeneralLinearMethod, T::Trees) =
    conjugate_order(underlying_one_step(glm, T).E, exact_solution(T), T)
