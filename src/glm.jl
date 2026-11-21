# Chapter 6: conjugacy of B-series, order by order.
#
# Two central B-series a and b (a(∅) = b(∅) = 1) are conjugate to order p
# if there is x with xa = bx up to trees of order p, i.e.
# x a x⁻¹ = b + O_{p+1}.  Algorithm 7 of the book takes a and b that
# already agree for |t| < p and decides whether one more order can be
# reached by a transformation x that differs from the identity only on
# trees of order p − 1.
#
# Book: Section 6.6, Algorithm 7 ("Determine whether two central B-series
# for which τ → 0 are conjugate").

"""
    conjugacy_step(a, b, T, p) -> (true, x a x⁻¹, x) or (false,)

Decide whether central B-series a and b with a = b + O_p are conjugate to
order p, i.e. whether there is x, differing from the identity only on the
trees of order p − 1, with x a x⁻¹ = b + O_{p+1}.

The values x(t), |t| = p − 1, are found from the conditions
(xa)(t ∘ τ) = (bx)(t ∘ τ).  These form a linear system, which is coupled
whenever a tree t ∘ τ can lose a leaf in more than one way; it is solved
exactly.  Then the conditions for all trees of order p are tested.

Book: Section 6.6, Algorithm 7 (corrected for the 2nd edition: the
1st edition computes each x(t) from t ∘ τ alone, which fails from p = 4).
"""
function conjugacy_step(a::AbstractVector, b::AbstractVector, T::Trees, p::Integer)
    unknowns = collect(trees_of_order(T, p - 1))
    rows = [beta_product(T, t, 1) for t in unknowns]
    residual(x) = compose(x, a, T) .- compose(b, x, T)
    r0 = residual(zero(a))
    M = zeros(eltype(a), length(rows), length(unknowns))
    for (j, t) in enumerate(unknowns)
        e = zero(a)
        e[t] = 1
        M[:, j] = residual(e)[rows] .- r0[rows]
    end
    sol = solve_exact(M, -r0[rows])
    sol === nothing && return (false,)
    x = zero(a)
    x[unknowns] = sol
    xa = compose(x, a, T)
    bx = compose(b, x, T)
    for t in trees_of_order(T, p)
        xa[t] == bx[t] || return (false,)
    end
    return (true, compose(xa, inverse(x, T), T), x)
end

# Solve M y = v by Gaussian elimination in exact arithmetic; free unknowns
# are set to 0 (for p = 2 the only unknown is x(τ), which is irrelevant,
# Book: Section 6.6, "The role of x₁").  Returns `nothing` if inconsistent.
function solve_exact(M::AbstractMatrix, v::AbstractVector; atol = 0)
    A = [M v]
    m, n = size(M)
    pivots = Int[]
    row = 1
    for col in 1:n
        row > m && break
        _, k = findmax(i -> abs(A[i, col]), row:m)   # partial pivoting
        k += row - 1
        abs(A[k, col]) > atol || continue
        A[[row, k], :] = A[[k, row], :]
        A[row, :] ./= A[row, col]
        for i in 1:m
            i != row && !iszero(A[i, col]) &&
                (A[i, :] .-= A[i, col] .* A[row, :])
        end
        push!(pivots, col)
        row += 1
    end
    all(abs(A[i, end]) ≤ atol for i in row:m) || return nothing
    y = zeros(eltype(A), n)
    for (i, col) in enumerate(pivots)
        y[col] = A[i, end]
    end
    return y
end
