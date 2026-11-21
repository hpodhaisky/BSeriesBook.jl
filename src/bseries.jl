# Chapter 3: B-series and their composition.
#
# A B-series a is stored as a vector over the standard numbering,
# a[n] = a(t_n), together with a(∅), which is 1 for the B-series of a
# numerical method (the "central" case).  The B-series itself is
#     B_h(a)y = a(∅) y + Σ_t h^|t| a(t)/σ(t) F(t)(y).
#
# Book: Sections 3.4–3.9; composition and Algorithms 5 and 6 in
# Section 3.9.

"""
    exact_solution(T) -> Vector{Rational{BigInt}}

The B-series E of the exact flow, E(t) = 1/t!.  Book: Section 3.5.
"""
exact_solution(T::Trees) = [1 // big(γ(T, n)) for n in 1:length(T)]

"""
    lambda_matrix(T, a, p_top) -> Matrix

The coefficients of λ(a, t) = Σ_{t' ≤ t} a(t ∖ t') t' (sum over the
subtrees t' of t that contain the root) for all trees of order ≤ `p_top`:
`λ[n, k]` is the coefficient of t_k in λ(a, t_n).

Book: Section 3.9, Algorithm 5 ("Evaluate λ"), using
λ(a, t₁ ∘ t₂) = a(t₂) λ(a, t₁) + λ(a, t₁) ∘ λ(a, t₂).
"""
function lambda_matrix(T::Trees, a::AbstractVector, p_top::Integer)
    (; last, L, R, prod) = T.table
    N = last[p_top]
    # [book:BSco10A:begin]
    λ = zeros(eltype(a), N, N)
    λ[1, 1] = 1
    for n in 2:N
        i = L[n]
        j = R[n]
        for k in 1:i
            λ[n, k] = a[j] * λ[i, k]
        end
        for k in 1:i, ℓ in 1:j
            if λ[i, k] != 0 && λ[j, ℓ] != 0
                λ[n, prod[k, ℓ]] = λ[n, prod[k, ℓ]] + λ[i, k] * λ[j, ℓ]
            end
        end
    end
    # [book:BSco10A:end]
    return λ
end

"""
    lambda_hat(T, a, p_top) -> Vector{Vector}

The efficient representation λ̂ of λ: `λ̂[n]` holds the coefficients of
λ(a, t_n) for the trees of order less than |t_n| (the coefficient of t_n
itself is 1, and no other tree of the same order occurs).

Book: Section 3.9, Algorithm 6 ("Evaluate λ", efficient form) and the
list printed after it.
"""
function lambda_hat(T::Trees, a::AbstractVector, p_top::Integer)
    (; order, last, L, R, prod) = T.table
    N = last[p_top]
    len(n) = order[n] == 1 ? 0 : last[order[n]-1]
    coeff(v, n, k) = k == n ? one(eltype(a)) : (k ≤ length(v[n]) ? v[n][k] : zero(eltype(a)))
    λ̂ = Vector{Vector{eltype(a)}}(undef, N)
    λ̂[1] = eltype(a)[]
    for n in 2:N
        i = L[n]
        j = R[n]
        λ̂[n] = zeros(eltype(a), len(n))
        # a(t_j) λ(a, t_i), including the term a(t_j) t_i
        for k in 1:i
            λ̂[n][k] += a[j] * coeff(λ̂, i, k)
        end
        # λ(a, t_i) ∘ λ(a, t_j) without the leading term t_i ∘ t_j = t_n
        for k in 1:i, ℓ in 1:j
            (k, ℓ) == (i, j) && continue
            c = coeff(λ̂, i, k) * coeff(λ̂, j, ℓ)
            c == 0 && continue
            λ̂[n][prod[k, ℓ]] += c
        end
    end
    return λ̂
end

"""
    compose(a, b, T; b0 = 1) -> Vector

The product ab of two B-series (first a, then b):
    (ab)(t) = b(∅) a(t) + Σ_{∅ < t' ≤ t} a(t ∖ t') b(t'),
for a with a(∅) = 1, where a(F) is the product of a over the trees of
the forest F.  Uses the table of [`splittings`](@ref).

Book: Section 3.9 (composition of B-series).
"""
function compose(a::AbstractVector, b::AbstractVector, T::Trees; b0 = 1)
    S = splittings(T)
    ab = b0 .* a
    for n in eachindex(S), (k, F, m) in S[n]
        ab[n] += m * b[k] * forest_value(a, F)
    end
    return ab
end

forest_value(a, F) = prod((a[f] for f in F); init = one(eltype(a)))

"""
    inverse(a, T) -> Vector

The inverse a⁻¹ in the B-group, from
    a⁻¹(t) = -a(t) - Σ_{∅ < t' < t} a(t ∖ t') a⁻¹(t').

Book: Section 3.9 ("Inverse of B-series").
"""
function inverse(a::AbstractVector, T::Trees)
    S = splittings(T)
    b = similar(a)
    for n in eachindex(S)
        b[n] = -a[n]
        for (k, F, m) in S[n]
            k == n && continue
            b[n] -= m * forest_value(a, F) * b[k]
        end
    end
    return b
end

"""
    identity_series(T) -> Vector

The identity of the B-group: id(∅) = 1, id(t) = 0.
"""
identity_series(T::Trees, ::Type{S} = Rational{BigInt}) where {S} =
    zeros(S, length(T))
