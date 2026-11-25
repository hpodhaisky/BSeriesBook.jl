# Chapter 3: B-series and their composition.
#
# A B-series is a mapping a: T# → ℝ, T# = {∅} ∪ T; its series is
#     B_h(a)y = a(∅) y + Σ_t h^|t| a(t)/σ(t) F(t)(y).
# The book (Section 3.4) uses three sets of such mappings:
#     𝖡*  all of them (a vector space),          GeneralBSeries,
#     𝖡   those with a(∅) = 1 (the B-group),      CentralBSeries,
#     𝖡⁰  those with a(∅) = 0 (a linear space),   VectorFieldBSeries.
# All three store a[n] = a(t_n) over the standard numbering together with
# their `Trees`; only GeneralBSeries stores a(∅).  The type is decided by
# the structure of a computation, never by the value of a(∅).
#
# Book: Sections 3.4–3.9; composition and Algorithms 5 and 6 in
# Section 3.9.

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    AbstractBSeries{S} <: AbstractVector{S}

A B-series a: T# → S, truncated at the order of its [`Trees`](@ref).  As a
vector, `a[n] = a(tₙ)` over the standard numbering; the value on the empty
tree is `a[∅]`.  Subtypes: [`GeneralBSeries`](@ref) (𝖡*),
[`CentralBSeries`](@ref) (𝖡, a(∅) = 1) and [`VectorFieldBSeries`](@ref)
(𝖡⁰, a(∅) = 0).

Indexing with a range, broadcasting and other generic vector operations
act on the coefficients a(tₙ) and return plain vectors; `+`, `-`,
scalar multiplication, `==` and `≈` take a(∅) into account and return
B-series.

Book: Section 3.4.
"""
abstract type AbstractBSeries{S} <: AbstractVector{S} end

"""
    GeneralBSeries(T, a0, coef)
    GeneralBSeries(f, T, a0)
    GeneralBSeries(a::AbstractBSeries)

A B-series in 𝖡* with a(∅) = `a0` and a(tₙ) = `coef[n]` (or a(t) = `f(t)`
for every [`Tree`](@ref) t).
Book: Section 3.4.
"""
struct GeneralBSeries{S} <: AbstractBSeries{S}
    T::Trees
    a0::S
    coef::Vector{S}
    function GeneralBSeries{S}(T::Trees, a0, coef::AbstractVector) where {S}
        check_length(T, coef)
        return new{S}(T, a0, coef)
    end
end

"""
    CentralBSeries(T, coef)
    CentralBSeries(f, T)
    CentralBSeries(a::AbstractBSeries)

A central B-series, an element of the B-group 𝖡: a(∅) = 1 (not stored)
and a(tₙ) = `coef[n]` (or a(t) = `f(t)` for every [`Tree`](@ref) t).  The
B-series of a one-step method.
Converting another B-series checks that a(∅) = 1;
`CentralBSeries{S}(a)` changes the number type to `S`.

Book: Section 3.4 ("central series").
"""
struct CentralBSeries{S} <: AbstractBSeries{S}
    T::Trees
    coef::Vector{S}
    function CentralBSeries{S}(T::Trees, coef::AbstractVector) where {S}
        check_length(T, coef)
        return new{S}(T, coef)
    end
end

"""
    VectorFieldBSeries(T, coef)
    VectorFieldBSeries(f, T)
    VectorFieldBSeries(a::AbstractBSeries)

A B-series in 𝖡⁰: a(∅) = 0 (not stored) and a(tₙ) = `coef[n]` (or
a(t) = `f(t)` for every [`Tree`](@ref) t).  `h⁻¹ B_h(a)` is a vector field, for example a modified
equation; the difference of two central B-series lies in 𝖡⁰.
Converting another B-series checks that a(∅) = 0.

Book: Section 3.4.
"""
struct VectorFieldBSeries{S} <: AbstractBSeries{S}
    T::Trees
    coef::Vector{S}
    function VectorFieldBSeries{S}(T::Trees, coef::AbstractVector) where {S}
        check_length(T, coef)
        return new{S}(T, coef)
    end
end

check_length(T, coef) = length(coef) == length(T) || throw(DimensionMismatch(
    "a B-series over $T needs $(length(T)) coefficients, got $(length(coef))"))

GeneralBSeries(T::Trees, a0, coef::AbstractVector) =
    GeneralBSeries{promote_type(typeof(a0), eltype(coef))}(T, a0, coef)
CentralBSeries(T::Trees, coef::AbstractVector) =
    CentralBSeries{eltype(coef)}(T, coef)
VectorFieldBSeries(T::Trees, coef::AbstractVector) =
    VectorFieldBSeries{eltype(coef)}(T, coef)

# From a function of the tree: f(t) for every t in T.
GeneralBSeries(f, T::Trees, a0) = GeneralBSeries(T, a0, [f(t) for t in T])
CentralBSeries(f, T::Trees) = CentralBSeries(T, [f(t) for t in T])
VectorFieldBSeries(f, T::Trees) = VectorFieldBSeries(T, [f(t) for t in T])

GeneralBSeries(a::AbstractBSeries) = GeneralBSeries(a.T, a[empty_tree], copy(a.coef))

# Change of the number type, e.g. CentralBSeries{BigFloat}(exact_solution(T)).
GeneralBSeries{S}(a::AbstractBSeries) where {S} =
    GeneralBSeries{S}(a.T, a[empty_tree], S.(a.coef))
CentralBSeries{S}(a::CentralBSeries) where {S} = CentralBSeries{S}(a.T, S.(a.coef))
VectorFieldBSeries{S}(a::VectorFieldBSeries) where {S} =
    VectorFieldBSeries{S}(a.T, S.(a.coef))
CentralBSeries(a::CentralBSeries) = copy(a)
VectorFieldBSeries(a::VectorFieldBSeries) = copy(a)
function CentralBSeries(a::AbstractBSeries)
    isone(a[empty_tree]) || throw(ArgumentError("a(∅) = $(a[∅]) ≠ 1: not a central B-series"))
    return CentralBSeries(a.T, copy(a.coef))
end
function VectorFieldBSeries(a::AbstractBSeries)
    iszero(a[empty_tree]) || throw(ArgumentError("a(∅) = $(a[∅]) ≠ 0: not in 𝖡⁰"))
    return VectorFieldBSeries(a.T, copy(a.coef))
end

# ---------------------------------------------------------------------------
# The empty tree and access
# ---------------------------------------------------------------------------

struct EmptyTree end

"""
    empty_tree

The empty tree ∅: `a[empty_tree]` is a(∅) for a B-series `a`.
Alias: `∅`.
"""
const empty_tree = EmptyTree()

const ∅ = empty_tree        # alias
Base.show(io::IO, ::EmptyTree) = print(io, "∅")

Base.size(a::AbstractBSeries) = size(a.coef)
Base.IndexStyle(::Type{<:AbstractBSeries}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(a::AbstractBSeries, n::Int) = a.coef[n]
Base.@propagate_inbounds Base.setindex!(a::AbstractBSeries, v, n::Int) =
    (a.coef[n] = v; a)
Base.getindex(a::AbstractBSeries, t::Tree) = a.coef[t.n]
Base.setindex!(a::AbstractBSeries, v, t::Tree) = (a.coef[t.n] = v; a)
Base.getindex(a::GeneralBSeries, ::EmptyTree) = a.a0
Base.getindex(::CentralBSeries{S}, ::EmptyTree) where {S} = one(S)
Base.getindex(::VectorFieldBSeries{S}, ::EmptyTree) where {S} = zero(S)

"""
    coefficients(a) -> Vector

The vector of the values a(tₙ), n = 1, …, length(T) (without a(∅)).  This
is the storage of `a`, not a copy.
"""
coefficients(a::AbstractBSeries) = a.coef

"""
    Trees(a::AbstractBSeries)

The trees over which the B-series `a` is defined.
"""
Trees(a::AbstractBSeries) = a.T

Base.copy(a::GeneralBSeries) = GeneralBSeries(a.T, a.a0, copy(a.coef))
Base.copy(a::CentralBSeries) = CentralBSeries(a.T, copy(a.coef))
Base.copy(a::VectorFieldBSeries) = VectorFieldBSeries(a.T, copy(a.coef))

# Two B-series can be combined if they are truncated at the same order;
# the standard numbering is then the same.
same_trees(a::AbstractBSeries, b::AbstractBSeries) =
    a.T.table.p_max == b.T.table.p_max
check_trees(a, b) = same_trees(a, b) || throw(DimensionMismatch(
    "B-series over different trees: $(a.T) and $(b.T)"))

Base.:(==)(a::AbstractBSeries, b::AbstractBSeries) =
    same_trees(a, b) && a[empty_tree] == b[empty_tree] && a.coef == b.coef
Base.isequal(a::AbstractBSeries, b::AbstractBSeries) =
    same_trees(a, b) && isequal(a[empty_tree], b[empty_tree]) && isequal(a.coef, b.coef)
Base.isapprox(a::AbstractBSeries, b::AbstractBSeries; kw...) =
    same_trees(a, b) && isapprox(a[empty_tree], b[empty_tree]; kw...) &&
    isapprox(a.coef, b.coef; kw...)

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", a::AbstractBSeries)
    T = a.T
    N = length(a)
    println(io, typeof(a), " to order ", T.table.p_max, " (", N, " trees)")
    rows = get(io, :limit, false) ? min(N, max(displaysize(io)[1] - 6, 5)) : N
    names = [bracket_string(T, n) for n in 1:rows]
    w = maximum(textwidth, names; init = 1)
    d = ndigits(rows)
    print(io, " ", " "^d, "  ", rpad("∅", w), "  ", a[empty_tree])
    for n in 1:rows
        print(io, "\n ", lpad(n, d), "  ", rpad(names[n], w), "  ", a[n])
    end
    rows < N && print(io, "\n  ⋮")
end

# ---------------------------------------------------------------------------
# Linear operations: the type follows from the type of the operands.
# A value of a(∅) that is fixed by the type is encoded as Val(0) or Val(1).
# ---------------------------------------------------------------------------

empty_value(a::GeneralBSeries) = a.a0
empty_value(::CentralBSeries) = Val(1)
empty_value(::VectorFieldBSeries) = Val(0)

fixed(x) = x == 0 ? Val(0) : x == 1 ? Val(1) : x
value(::Val{x}) where {x} = x
value(x) = x
add_empty(x, y) = value(x) + value(y)
add_empty(::Val{x}, ::Val{y}) where {x,y} = fixed(x + y)
sub_empty(x, y) = value(x) - value(y)
sub_empty(::Val{x}, ::Val{y}) where {x,y} = fixed(x - y)
mul_empty(c, x) = c * value(x)
mul_empty(c, ::Val{0}) = Val(0)
div_empty(x, c, S) = S(value(x)) / c
div_empty(::Val{0}, c, S) = Val(0)

make(T, ::Val{0}, coef) = VectorFieldBSeries(T, coef)
make(T, ::Val{1}, coef) = CentralBSeries(T, coef)
make(T, a0, coef) = GeneralBSeries(T, a0, coef)

function Base.:+(a::AbstractBSeries, b::AbstractBSeries)
    check_trees(a, b)
    return make(a.T, add_empty(empty_value(a), empty_value(b)), a.coef + b.coef)
end
function Base.:-(a::AbstractBSeries, b::AbstractBSeries)
    check_trees(a, b)
    return make(a.T, sub_empty(empty_value(a), empty_value(b)), a.coef - b.coef)
end
Base.:-(a::AbstractBSeries) = make(a.T, mul_empty(-1, empty_value(a)), -a.coef)
Base.:*(c::Number, a::AbstractBSeries) = make(a.T, mul_empty(c, empty_value(a)), c * a.coef)
Base.:*(a::AbstractBSeries, c::Number) = c * a
Base.:/(a::AbstractBSeries, c::Number) = make(a.T, div_empty(empty_value(a), c, eltype(a)), a.coef / c)

"""
    zero(a::AbstractBSeries) -> VectorFieldBSeries

The zero element of 𝖡* (it lies in 𝖡⁰).
"""
Base.zero(a::AbstractBSeries) = VectorFieldBSeries(a.T, zero(a.coef))

# ---------------------------------------------------------------------------
# Special members
# ---------------------------------------------------------------------------

"""
    exact_solution(T) -> CentralBSeries{Rational{BigInt}}

The B-series E of the exact flow, E(t) = 1/t!.  Book: Section 3.5.
"""
exact_solution(T::Trees) = CentralBSeries(t -> 1 // big(factorial(t)), T)

"""
    identity_series(T, S = Rational{BigInt}) -> CentralBSeries{S}

The identity of the B-group: id(∅) = 1, id(t) = 0.  Also `one(a)`.
"""
identity_series(T::Trees, ::Type{S} = Rational{BigInt}) where {S} =
    CentralBSeries(T, zeros(S, length(T)))

Base.one(a::CentralBSeries) = identity_series(a.T, eltype(a))

"""
    rescale(a, θ) -> B-series of the same type

The B-series of the step size θh: a(t) θ^|t|.  For the exact solution,
`rescale(E, θ)` is the exact flow over θh.
"""
rescale(a::AbstractBSeries, theta) =
    make(a.T, empty_value(a), [a[n] * theta^a.T.table.order[n] for n in 1:length(a)])

# ---------------------------------------------------------------------------
# λ (Algorithms 5 and 6)
# ---------------------------------------------------------------------------

"""
    lambda_matrix(a::CentralBSeries, p_top = p_max) -> Matrix

The coefficients of λ(a, t) = Σ_{t' ≤ t} a(t ∖ t') t' (sum over the
subtrees t' of t that contain the root) for all trees of order ≤ `p_top`:
`λ[n, k]` is the coefficient of t_k in λ(a, t_n), so that
(ab)(tₙ) = b(∅) a(tₙ) + Σ_k λ[n, k] b(t_k).

Book: Section 3.9, Algorithm 5 ("Evaluate λ"), using
λ(a, t₁ ∘ t₂) = a(t₂) λ(a, t₁) + λ(a, t₁) ∘ λ(a, t₂).
"""
function lambda_matrix(a::CentralBSeries, p_top::Integer = a.T.table.p_max)
    (; last, L, R, prod) = a.T.table
    N = last[p_top]
    # [book:BSco10A:begin]
    lambda = zeros(eltype(a), N, N)
    lambda[1, 1] = 1
    for n in 2:N
        i = L[n]
        j = R[n]
        for k in 1:i
            lambda[n, k] = a[j] * lambda[i, k]
        end
        for k in 1:i, l in 1:j
            if lambda[i, k] != 0 && lambda[j, l] != 0
                lambda[n, prod[k, l]] = lambda[n, prod[k, l]] + lambda[i, k] * lambda[j, l]
            end
        end
    end
    # [book:BSco10A:end]
    return lambda
end

"""
    lambda_hat(a::CentralBSeries, p_top = p_max) -> Vector{Vector}

The efficient representation λ̂ of λ: `λ̂[n]` holds the coefficients of
λ(a, t_n) for the trees of order less than |t_n| (the coefficient of t_n
itself is 1, and no other tree of the same order occurs).

Book: Section 3.9, Algorithm 6 ("Evaluate λ", efficient form) and the
list printed after it.
"""
function lambda_hat(a::CentralBSeries, p_top::Integer = a.T.table.p_max)
    (; order, last, L, R, prod) = a.T.table
    S = eltype(a)
    N = last[p_top]
    len(n) = order[n] == 1 ? 0 : last[order[n]-1]
    coeff(v, n, k) = k == n ? one(S) : (k <= length(v[n]) ? v[n][k] : zero(S))
    lhat = Vector{Vector{S}}(undef, N)
    lhat[1] = S[]
    for n in 2:N
        i = L[n]
        j = R[n]
        lhat[n] = zeros(S, len(n))
        # a(t_j) λ(a, t_i), including the term a(t_j) t_i
        for k in 1:i
            lhat[n][k] += a[j] * coeff(lhat, i, k)
        end
        # λ(a, t_i) ∘ λ(a, t_j) without the leading term t_i ∘ t_j = t_n
        for k in 1:i, l in 1:j
            (k, l) == (i, j) && continue
            c = coeff(lhat, i, k) * coeff(lhat, j, l)
            c == 0 && continue
            lhat[n][prod[k, l]] += c
        end
    end
    return lhat
end

# ---------------------------------------------------------------------------
# Composition and inverse
# ---------------------------------------------------------------------------

"""
    compose(a::CentralBSeries, b::AbstractBSeries), a * b

The product ab of the book (first a, then b), 𝖡 × 𝖡* → 𝖡*:
    (ab)(∅) = b(∅),
    (ab)(t) = b(∅) a(t) + Σ_{∅ < t' ≤ t} a(t ∖ t') b(t'),
where a(F) is the product of a over the trees of the forest F.  The
result has the type of `b`: it is central if b is central and in 𝖡⁰ if b
is.  The first factor must be central, otherwise B_h(b) ∘ B_h(a) is not a
B-series.  Uses the table of [`splittings`](@ref).

Book: Section 3.9 (composition of B-series).
"""
function compose(a::CentralBSeries, b::AbstractBSeries)
    check_trees(a, b)
    return make(a.T, empty_value(b),
                compose_coefficients(a.coef, b.coef, b[empty_tree], a.T))
end

Base.:*(a::CentralBSeries, b::AbstractBSeries) = compose(a, b)
Base.:^(a::CentralBSeries, p::Integer) =
    p < 0 ? Base.power_by_squaring(inverse(a), -p) : Base.power_by_squaring(a, p)

# Kernel on coefficient vectors: a(∅) = 1, b(∅) = b0.
function compose_coefficients(a::AbstractVector, b::AbstractVector, b0, T::Trees)
    S = splitting_table(T)
    ab = promote_type(eltype(a), eltype(b), typeof(b0)).(b0 .* a)
    for n in eachindex(S), (k, F, m) in S[n]
        ab[n] += m * b[k] * forest_value(a, F)
    end
    return ab
end

forest_value(a, F) = prod((a[f] for f in F); init = one(eltype(a)))

"""
    slope(y::CentralBSeries, gamma = 0) -> VectorFieldBSeries

The B-series of hf(w), where w = y + γ hf(w) (as B-series about y₀):
returns dw = wD with w = y + γ dw, from
    dw(τ) = 1,   dw(tₙ) = dw(t_{L[n]}) w(t_{R[n]}),   w(tₙ) = y(tₙ) + γ dw(tₙ),
in one loop over the standard numbering, without iteration.

`slope(y)` (γ = 0) is yD, the series of hf(B_h(y)y₀);
`slope(identity_series(T))` is D itself.  With y = id and γ = 1,
w is the implicit Euler method.

Book: Section 3.5 (the map slope, D ∈ 𝖡⁰, Theorem on aD).
"""
function slope(y::CentralBSeries, gamma = 0)
    (; L, R) = y.T.table
    S = promote_type(eltype(y), typeof(gamma))
    dw = zeros(S, length(y))
    w = zeros(S, length(y))
    for n in eachindex(dw)
        dw[n] = n == 1 ? one(S) : dw[L[n]] * w[R[n]]
        w[n] = y[n] + gamma * dw[n]
    end
    return VectorFieldBSeries(y.T, dw)
end

"""
    inverse(a::CentralBSeries), inv(a)

The inverse a⁻¹ in the B-group, from
    a⁻¹(t) = -a(t) - Σ_{∅ < t' < t} a(t ∖ t') a⁻¹(t').

Book: Section 3.9 ("Inverse of B-series").
"""
inverse(a::CentralBSeries) = CentralBSeries(a.T, inverse_coefficients(a.coef, a.T))

Base.inv(a::CentralBSeries) = inverse(a)

function inverse_coefficients(a::AbstractVector, T::Trees)
    S = splitting_table(T)
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
    order_of_accuracy(a::CentralBSeries; atol = 0) -> Int

The largest p ≤ p_max with a(t) = E(t) = 1/t! for all trees with |t| ≤ p.
Use `atol > 0` for floating-point coefficients.

Book: Section 3.6 (order conditions).
"""
function order_of_accuracy(a::CentralBSeries; atol = 0)
    T = a.T
    for p in 1:T.table.p_max
        for n in tree_range(T, p)
            isapprox(a[n], 1 // T.fns.factorial[n]; atol = atol) || return p - 1
        end
    end
    return T.table.p_max
end
