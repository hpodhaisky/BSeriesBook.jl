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
#   5. the modified equation of the classical RK method (rationals);
# at the end the one-off costs (tables of splittings and partitions) and
# further operations of BSeriesBook.jl alone.
# BSeriesBook times exclude building the tree table `Trees(p)` (reported
# separately in line 1) and the table of splittings, which is built once on
# the first call of `compose`, and the table of partitions, built once on
# the first call of `modified_equation` (both listed at the end).

using BSeriesBook, BenchmarkTools
import BSeries, RootedTrees as RT

rt(t) = RT.rootedtree(level_sequence(t))
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
    trees = [rt(t) for t in T]
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
    println("trees+sigma+t!   Book ", fmt(tm(t1)), "   RootedTrees ", fmt(tm(t1b)))
    # 2. RK B-series (rational)
    Phi = elementary_weights(rk, T); bs = BSeries.bseries(rkBS, p)
    @assert all(bs[trees[n]] == Phi[n] for n in 1:N)
    t2 = @benchmark elementary_weights($rk, $T) samples=20 evals=1
    t2b = @benchmark BSeries.bseries($rkBS, $p) samples=20 evals=1
    println("RK B-series (Q)  Book ", fmt(tm(t2)), "   BSeries ", fmt(tm(t2b)))
    t2f = @benchmark elementary_weights($rkF, $T) samples=20 evals=1
    t2fb = @benchmark BSeries.bseries($rkBSF, $p) samples=20 evals=1
    println("RK B-series (F64)Book ", fmt(tm(t2f)), "   BSeries ", fmt(tm(t2fb)))
    # 3. composition (Float64)
    a = CentralBSeries(T, rand(N)); b = CentralBSeries(T, rand(N))
    A = toBS(a, trees); B = toBS(b, trees)
    ab = compose(a, b); AB = BSeries.compose(A, B)
    @assert maximum(abs(AB[trees[n]] - ab[n]) for n in 1:N) < 1e-9 * maximum(abs, ab)
    t3 = @benchmark compose($a, $b) samples=5 evals=1
    t3b = @benchmark BSeries.compose($A, $B) samples=5 evals=1
    println("composition (F64)Book ", fmt(tm(t3)), "   BSeries ", fmt(tm(t3b)))
    # 4. order of accuracy
    t4 = @benchmark order_of_accuracy($rk, $T) samples=20 evals=1
    t4b = @benchmark BSeries.order_of_accuracy(BSeries.bseries($rkBS, $p)) samples=20 evals=1
    println("order (Q)        Book ", fmt(tm(t4)), "   BSeries ", fmt(tm(t4b)))
    # 5. modified equation (rational)
    Phib = CentralBSeries(T, Rational{BigInt}.(Phi)); bsb = BSeries.bseries(BSeries.RungeKuttaMethod(big.(rk.A), big.(rk.b)), p)
    me = modified_equation(Phib); meb = BSeries.modified_equation(bsb)
    @assert all(meb[trees[n]] == me[n] for n in 1:N)
    t5 = @benchmark modified_equation($Phib) samples=5 evals=1
    t5b = @benchmark BSeries.modified_equation($bsb) samples=3 evals=1
    println("modified eq. (Q) Book ", fmt(tm(t5)), "   BSeries ", fmt(tm(t5b)))
end

# ---------------------------------------------------------------------------
# One-off costs and further operations (BSeriesBook.jl only); the GLM is
# method (6.5 h)(ii) of the book.  Each table is built on a fresh Trees(p).
# ---------------------------------------------------------------------------
glm = GeneralLinearMethod(
    [0 0 0 0; 3//4 0 0 0; 93//250 -9//125 0 0; 8881//8640 1069//3240 -1855//5184 0],
    [1 0; 1 1; 1 44//125; 1 691//540],
    [5//54 32//81 250//567 1//14; -19//16 -37//6 1175//336 27//7],
    [1 0; 0 0])
println("\n=== one-off costs and further operations")
for p in (6, 8, 10)
    T = Trees(p); N = length(T)
    nsplit = sum(length, splittings(T))
    tsplit = @belapsed splittings(T0) setup = (T0 = Trees($p)) evals = 1 samples = 10
    a = CentralBSeries(T, rand(N))
    tlam = @belapsed lambda_matrix($a)
    tinv = @belapsed inverse($a)
    println("p = $p: trees $N, splittings $nsplit, building splittings ", fmt(tsplit),
            ", λ matrix (F64) ", fmt(tlam), ", inverse (F64) ", fmt(tinv))
    p == 6 && continue
    npart = sum(length, tree_partitions(T))
    tpart = @belapsed tree_partitions(T0) setup = (T0 = Trees($p)) evals = 1 samples = 3
    tglm = @belapsed underlying_one_step($glm, $T) samples = 5 evals = 1
    tord = @belapsed order_of_accuracy($glm, $T) samples = 3 evals = 1
    println("        partitions $npart, building partitions ", fmt(tpart),
            ", underlying one-step method (Q) ", fmt(tglm),
            ", order_of_accuracy(glm) ", fmt(tord))
end
