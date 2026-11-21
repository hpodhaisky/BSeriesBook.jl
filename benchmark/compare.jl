# Compare BSeriesBook with BSeries.jl / RootedTrees.jl: results and run time.
#
# Run from this directory:
#     julia --project=. compare.jl
#
# For each order p = 6, 8, 10 the script first checks that both packages
# give the same results, then reports the minimum run time of
#   1. generating all trees with σ and t!,
#   2. the B-series of the classical RK method (rationals and Float64),
#   3. the composition of two B-series (Float64),
#   4. the order of the classical RK method,
#   5. the modified equation of the classical RK method (rationals).
# BSeriesBook times exclude building the tree table `Trees(p)` (reported
# separately in line 1) and the table of splittings, which is built once on
# the first call of `compose` (11 ms at order 10), and the table of
# partitions, built once on the first call of `modified_equation`
# (0.1 s at order 10).

using BSeriesBook, BenchmarkTools
import BSeries, RootedTrees as RT

rt(T, n) = RT.rootedtree(level_sequence(T, n))
function toBS(v, trees)
    s = BSeries.TruncatedBSeries{eltype(trees),eltype(v)}()
    s[RT.rootedtree(Int[])] = one(eltype(v))
    for n in eachindex(trees); s[trees[n]] = v[n]; end
    s
end
rk = kutta4(); rkBS = BSeries.RungeKuttaMethod(rk.A, rk.b)
rkF = RungeKuttaMethod(Float64.(rk.A), Float64.(rk.b)); rkBSF = BSeries.RungeKuttaMethod(Float64.(rk.A), Float64.(rk.b))
fmt(t) = t < 1e-3 ? string(round(t*1e6, digits=1), " µs") : t < 1 ? string(round(t*1e3, digits=2), " ms") : string(round(t, digits=2), " s")
tm(b) = minimum(b).time / 1e9

for p in (6, 8, 10)
    T = Trees(p); N = length(T)
    trees = [rt(T, n) for n in 1:N]
    println("\n=== order $p, $N trees")
    # 1. trees + tree functions
    t1 = @benchmark Trees($p) samples=20 evals=1
    t1b = @benchmark begin
        s = 0
        for q in 1:$p, t in RT.RootedTreeIterator(q)
            s += RT.symmetry(t) + RT.density(t)
        end
        s
    end samples=20 evals=1
    println("trees+σ+γ        Book ", fmt(tm(t1)), "   RootedTrees ", fmt(tm(t1b)))
    # 2. RK B-series (rational)
    Φ = elementary_weights(rk, T); bs = BSeries.bseries(rkBS, p)
    @assert all(bs[trees[n]] == Φ[n] for n in 1:N)
    t2 = @benchmark elementary_weights($rk, $T) samples=20 evals=1
    t2b = @benchmark BSeries.bseries($rkBS, $p) samples=20 evals=1
    println("RK B-series (Q)  Book ", fmt(tm(t2)), "   BSeries ", fmt(tm(t2b)))
    t2f = @benchmark elementary_weights($rkF, $T) samples=20 evals=1
    t2fb = @benchmark BSeries.bseries($rkBSF, $p) samples=20 evals=1
    println("RK B-series (F64)Book ", fmt(tm(t2f)), "   BSeries ", fmt(tm(t2fb)))
    # 3. composition (Float64)
    a = rand(N); b = rand(N)
    A = toBS(a, trees); B = toBS(b, trees)
    ab = compose(a, b, T); AB = BSeries.compose(A, B)
    @assert maximum(abs(AB[trees[n]] - ab[n]) for n in 1:N) < 1e-9 * maximum(abs, ab)
    t3 = @benchmark compose($a, $b, $T) samples=5 evals=1
    t3b = @benchmark BSeries.compose($A, $B) samples=5 evals=1
    println("composition (F64)Book ", fmt(tm(t3)), "   BSeries ", fmt(tm(t3b)))
    # 4. order of accuracy
    t4 = @benchmark order_of_accuracy($rk, $T) samples=20 evals=1
    t4b = @benchmark BSeries.order_of_accuracy(BSeries.bseries($rkBS, $p)) samples=20 evals=1
    println("order (Q)        Book ", fmt(tm(t4)), "   BSeries ", fmt(tm(t4b)))
    # 5. modified equation (rational)
    Φb = Rational{BigInt}.(Φ); bsb = BSeries.bseries(BSeries.RungeKuttaMethod(big.(rk.A), big.(rk.b)), p)
    me = modified_equation(Φb, T); meb = BSeries.modified_equation(bsb)
    @assert all(meb[trees[n]] == me[n] for n in 1:N)
    t5 = @benchmark modified_equation($Φb, $T) samples=5 evals=1
    t5b = @benchmark BSeries.modified_equation($bsb) samples=3 evals=1
    println("modified eq. (Q) Book ", fmt(tm(t5)), "   BSeries ", fmt(tm(t5b)))
end
