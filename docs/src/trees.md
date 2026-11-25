# Standard numbering

The book numbers the rooted trees by Algorithm 3 ("Generate trees",
Section 2.5). Every tree except
`t₁ = τ` is the β-product `tₙ = t_{L[n]} ∘ t_{R[n]}`: the tree `t_{R[n]}`
is grafted onto the root of `t_{L[n]}`, and `R[n] ≤ R[L[n]]` makes the
decomposition unique. Within each order the trees are numbered in the
order in which the algorithm produces them. In the package a tree is a
[`Tree`](@ref), `T[n]`, which carries its standard number `n`
([`number`](@ref)) but is not an integer.

The tables on this page are generated from the package when the
documentation is built ([`Trees`](@ref), [`generate_trees`](@ref)). The
drawings have the root at the bottom, as in the book.

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

## Trees up to order 5

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

## Product table

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

## In the REPL

```@example trees
T = Trees(5)
[(t, T.table.L[number(t)], T.table.R[number(t)]) for t in trees_of_order(T, 5)]
```

Larger tables are obtained the same way, e.g. `Trees(10)` has 1205 trees.
The level sequence of a tree ([`level_sequence`](@ref)) identifies it in
RootedTrees.jl.
