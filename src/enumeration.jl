# Chapter 2: counting trees.
#
# Book: Section 2.4, Algorithms 1 and 2, the table of tree enumerations,
# Theorems 2.4B and 2.4E.

"""
    count_trees(p_max) -> Vector{BigInt}

`Ntrees[p]`, the number of rooted trees of order `p = 1:p_max`.

Book: Section 2.4, Algorithm 1 ("Find the number of trees of each
order"), based on Theorem 2.4B.

```jldoctest
julia> count_trees(8)
8-element Vector{BigInt}:
   1
   1
   2
   4
   9
  20
  48
 115
```
"""
function count_trees(p_max::Integer)
    # [book:TRsp10A:begin]
    Ntrees = ones(BigInt, p_max)
    for i in 2:p_max-1
        for j in p_max:-1:1
            a = Ntrees[i]
            for k in 1:(j-1)÷i
                Ntrees[j] = Ntrees[j] + a * Ntrees[j-i*k]
                a = a * (Ntrees[i] + k) ÷ (k + 1)
            end
        end
    end
    # [book:TRsp10A:end]
    return Ntrees
end

"""
    count_unrooted_trees(p_max) -> (Utrees, UNtrees)

`Utrees[p]`, the number of unrooted (free) trees of order `p`, and
`UNtrees[p]`, the number of non-superfluous unrooted trees.

Book: Section 2.4, Algorithm 2 ("Find the number of unrooted and
non-superfluous unrooted trees of each order"), based on Theorem 2.4E.
"""
function count_unrooted_trees(p_max::Integer)
    Ntrees = count_trees(p_max)
    # [book:TRsp20A:begin]
    Xtrees = zeros(BigInt, p_max)
    Ytrees = zeros(BigInt, p_max)
    for i in 1:p_max-1
        for j in 1:p_max-i
            Xtrees[i+j] = Xtrees[i+j] + Ntrees[i] * Ntrees[j]
        end
    end
    for i in 1:p_max÷2
        Ytrees[2*i] = Ntrees[i]
    end
    Utrees = Ntrees - (Xtrees - Ytrees) .÷ 2
    UNtrees = Ntrees - (Xtrees + Ytrees) .÷ 2
    # [book:TRsp20A:end]
    return Utrees, UNtrees
end
