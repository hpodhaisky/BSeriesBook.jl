# Chapter 3: substitution of B-series and the modified equation.
#
# A partition of a tree t is a set of edges that is cut; it splits t into
# pieces (subtrees), and contracting every piece to a single vertex gives
# the skeleton.  For a B-series b of a vector field (b(∅) = 0) and a
# B-series a, the substitution law of Chartier, Hairer and Vilmart is
#     (b ⋆ a)(t) = Σ_partitions  b(piece₁) ⋯ b(piece_m) a(skeleton),
# with the same normalisation  Σ h^|t| a(t)/σ(t) F(t)  as elsewhere:
# b ⋆ a is the B-series obtained when the method with B-series a is
# applied to the (modified) vector field  h⁻¹ B_h(b).
#
# Book: Section 3.9 (substitution); P. Chartier, E. Hairer, G. Vilmart,
# Algebraic structures of B-series, Found. Comput. Math. 10 (2010).

"""
    join_roots(s::Tree, t::Tree) -> Tree

The tree obtained from s and t by identifying their roots,
`[B⁻(s) B⁻(t)]` (Butcher's product s × t).
"""
function join_roots(s::Tree, t::Tree)
    T = s.T.table.p_max >= t.T.table.p_max ? s.T : t.T
    order(s) + order(t) - 1 <= T.table.p_max || throw(DomainError(
        (s, t), "the tree has order > $(T.table.p_max)"))
    return Tree(T, join_root_numbers(T, s.n, t.n))
end

join_root_numbers(T::Trees, s, u) =
    tree_number(T, [child_numbers(T, s); child_numbers(T, u)])

"""
    tree_partitions(T) -> Vector{Vector{Partition}}

For every tree t_n all partitions `(k, F, s, m)`: root piece t_k, other
pieces F, skeleton t_s, multiplicity m.  Built once per `Trees` object by
the β-product recursion: for t_n = t_i ∘ t_j combine every partition
(k, F, s) of t_i with every partition (ℓ, G, u) of t_j;
- the edge to t_j is cut:  root piece t_k, pieces F ∪ G ∪ {t_ℓ},
  skeleton t_s ∘ t_u;
- it is kept:  root piece t_k ∘ t_ℓ, pieces F ∪ G, skeleton t_s × t_u
  (roots identified, [`join_roots`](@ref)).
Summed over F, the multiplicities of t_n add up to 2^(|t_n| - 1).
"""
function tree_partitions(T::Trees)
    T.parts[] === nothing || return T.parts[]
    (; L, R, prod) = T.table
    P = Vector{Vector{Partition}}(undef, length(T))
    P[1] = [(1, Int[], 1, 1)]
    for n in 2:length(T)
        i, j = L[n], R[n]
        mult = Dict{Tuple{Int,Vector{Int},Int},Int}()
        add!(k, F, s, m) =
            (key = (k, sort(F), s); mult[key] = get(mult, key, 0) + m)
        for (k, F, s, m) in P[i], (l, G, u, m2) in P[j]
            add!(k, [F; G; l], prod[s, u], m * m2)          # edge cut
            add!(prod[k, l], [F; G], join_root_numbers(T, s, u), m * m2)
        end
        P[n] = sort!([(k, F, s, m) for ((k, F, s), m) in mult])
    end
    T.parts[] = P
    return P
end

"""
    substitute(b::VectorFieldBSeries, a::AbstractBSeries)

The substitution b ⋆ a of the B-series b ∈ 𝖡⁰ of a vector field into the
B-series a: the method a applied to the vector field h⁻¹ B_h(b),
    (b ⋆ a)(∅) = a(∅),
    (b ⋆ a)(t) = Σ_partitions b(pieces) a(skeleton).
The result has the type of `a`.  Same argument order as
`BSeries.substitute`.

Book: Section 3.9; Chartier–Hairer–Vilmart (2010), Section 3.2.
"""
function substitute(b::VectorFieldBSeries, a::AbstractBSeries)
    check_trees(a, b)
    P = tree_partitions(a.T)
    c = zeros(promote_type(eltype(a), eltype(b)), length(P))
    for n in eachindex(P), (k, F, s, m) in P[n]
        c[n] += m * b[k] * forest_value(b, F) * a[s]
    end
    return make(a.T, empty_value(a), c)
end

"""
    modified_equation(a::CentralBSeries) -> VectorFieldBSeries

The B-series b of the modified vector field of a method with B-series a:
the exact flow of y′ = h⁻¹ B_h(b)(y) is the method, b ⋆ E = a.  Solved
tree by tree (in standard order) from
    b(t) = a(t) − Σ_{partitions other than {t}} b(pieces) E(skeleton),
the partition {t} (no edge cut, skeleton τ) being the only one in which
b(t) occurs.  For a method of order p, b(t) = 0 for 1 < |t| ≤ p.

Book: backward error analysis (not in the 1st edition); Hairer, Lubich,
Wanner, Geometric Numerical Integration, Chapter IX.
"""
modified_equation(a::CentralBSeries) = solve_substitution(a, exact_solution(a.T))

"""
    modified_equation(rk::RungeKuttaMethod, T)

The modified vector field of the Runge–Kutta method `rk`.
"""
modified_equation(rk::RungeKuttaMethod, T::Trees) =
    modified_equation(bseries(rk, T))

"""
    modifying_integrator(a::CentralBSeries) -> VectorFieldBSeries
    modifying_integrator(rk::RungeKuttaMethod, T)

The B-series b of the vector field to which the method with B-series a
must be applied to give the exact flow: b ⋆ a = E.  Solved tree by tree
from b(t) a(τ) = E(t) − Σ_{other partitions} b(pieces) a(skeleton);
requires a(τ) ≠ 0.
"""
modifying_integrator(a::CentralBSeries) = solve_substitution(exact_solution(a.T), a)

modifying_integrator(rk::RungeKuttaMethod, T::Trees) =
    modifying_integrator(bseries(rk, T))

# Solve b ⋆ e = target for b ∈ 𝖡⁰, tree by tree.  Every partition other
# than the trivial one (k = n, no cut) involves b only on trees of lower
# order.
function solve_substitution(target::CentralBSeries, e::CentralBSeries)
    check_trees(target, e)
    P = tree_partitions(e.T)
    S = promote_type(eltype(target), eltype(e))
    b = zeros(S, length(P))
    for n in eachindex(P)
        r = S(target[n])
        for (k, F, s, m) in P[n]
            k == n && continue
            r -= m * b[k] * forest_value(b, F) * e[s]
        end
        b[n] = r / e[1]
    end
    return VectorFieldBSeries(e.T, b)
end
