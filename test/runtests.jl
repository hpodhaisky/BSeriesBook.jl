using Test
using BSeriesBook
import BSeries
import RootedTrees as RT

include("book_values.jl")

# Map from the standard numbering to RootedTrees.jl trees (level sequences)
rt(t::Tree) = RT.rootedtree(level_sequence(t))
rt(T::Trees, n::Integer) = rt(T[n])        # entries of the splitting/partition tables

randq(rng_state) = (rng_state[] = (1103515245 * rng_state[] + 12345) % 2^31;
                    (rng_state[] % 19 - 9) // (rng_state[] % 7 + 1))
randc(T, state) = CentralBSeries(_ -> big(randq(state)), T)
randv(T, state) = VectorFieldBSeries(_ -> big(randq(state)), T)

@testset "BSeriesBook" begin

@testset "Ch. 2: Algorithm 3, generate trees (book values to order 5)" begin
    T = generate_trees(5)
    @test length(T) == 17
    @test T.L[2:17] == BOOK_L
    @test T.R[2:17] == BOOK_R
    for (l, row) in BOOK_prod, (r, n) in enumerate(row)
        @test T.prod[l, r] == n
    end
end

@testset "Ch. 2: Algorithm 4, tree functions (book values to order 5)" begin
    F = tree_functions(generate_trees(5))
    @test F.symmetry == BOOK_symmetry
    @test F.factorial == BOOK_factorial
    @test F.Bminus == BOOK_Bminus
end

@testset "Ch. 2: Algorithms 1 and 2, counting trees" begin
    @test count_trees(15) == BOOK_a
    U, UN = count_unrooted_trees(15)
    @test U == BOOK_b
    @test UN == BOOK_c
    @test count_trees(10) == [RT.count_trees(p) for p in 1:10]
end

@testset "Ch. 2: trees and tree functions vs. RootedTrees.jl (order 8)" begin
    T = Trees(8)
    @test length(T) == sum(count_trees(8))
    seen = Set{Any}()
    for n in 1:length(T)
        t = rt(T[n])
        push!(seen, t)
        @test RT.order(t) == order(T[n])
        @test RT.symmetry(t) == symmetry(T[n])
        @test RT.density(t) == factorial(T[n])
        @test RT.α(t) == alpha(T[n])       # RootedTrees.jl has only the Greek names
        @test RT.β(t) == beta(T[n])
        @test bplus(T, children(T[n])) == T[n]
        @test number(T[n]) == n
    end
    @test length(seen) == length(T)            # all trees distinct
    # β-product: t_ℓ ∘ t_r is the Butcher product `∘` in RootedTrees.jl
    for l in 1:length(T), r in 1:length(T)
        order(T[l]) + order(T[r]) <= 8 || continue
        @test rt(beta_product(T[l], T[r])) ==
              rt(T[l]) ∘ rt(T[r])          # ∘: Butcher product of RootedTrees.jl
    end
end

@testset "Ch. 2: the type Tree" begin
    T = Trees(5)
    t = T[17]
    @test t isa Tree && number(t) == 17 && order(t) == 5
    @test bracket(t) == "[[[[τ]]]]" && sprint(show, t) == "t17 = [[[[τ]]]]"
    @test factorial(T[8]) == 24 && factorial(4) == 24     # tree and integer
    @test collect(T) == [T[n] for n in 1:length(T)] && first(T) == T[1]
    @test trees_of_order(T, 3) == [T[3], T[4]]
    @test beta_product(T[2], T[1]) == T[3] && join_roots(T[2], T[2]) == T[3]
    @test bplus(T, Tree[]) == T[1] && children(T[1]) == Tree[]
    @test T[17] == Trees(8)[17]              # the numbering does not depend on p_max
    @test_throws BoundsError T[18]
    @test_throws BoundsError T[0]
    @test_throws DomainError beta_product(T[17], T[1])
    @test_throws DomainError bplus(T, [T[17]])
    @test_throws MethodError order(17)
    @test_throws MethodError order(T, 17)
    a = exact_solution(T)
    @test a[t] == 1 // 120 && a[T[1]] == 1
    a[t] = 7
    @test a[17] == 7
end

@testset "Ch. 3: the types 𝖡*, 𝖡, 𝖡⁰" begin
    T = Trees(5)
    state = Ref(2)
    a, b = randc(T, state), randc(T, state)
    v, w = randv(T, state), randv(T, state)
    g = GeneralBSeries(_ -> big(randq(state)), T, big(5))
    @test a[empty_tree] == 1 && v[empty_tree] == 0 && g[empty_tree] == 5 && a[empty_tree] == 1
    @test length(a) == length(T) && a[3] == coefficients(a)[3]
    @test Trees(a) === T
    # the type of a linear combination follows from a(∅), fixed by the types
    @test a - b isa VectorFieldBSeries
    @test a + v isa CentralBSeries && v + a isa CentralBSeries
    @test a - v isa CentralBSeries
    @test v + w isa VectorFieldBSeries && 3v isa VectorFieldBSeries && -v isa VectorFieldBSeries
    @test a + b isa GeneralBSeries && (a + b)[empty_tree] == 2
    @test v - a isa GeneralBSeries && (v - a)[empty_tree] == -1
    @test 2a isa GeneralBSeries && (a / 2)[empty_tree] == 1 // 2
    @test eltype(a / 2) == Rational{BigInt} && v / 2 isa VectorFieldBSeries
    @test g + v isa GeneralBSeries && (g - a)[empty_tree] == 4
    @test zero(a) isa VectorFieldBSeries && a + zero(a) == a
    # type stable
    @test @inferred(a - b) isa VectorFieldBSeries
    @test @inferred(a + b) isa GeneralBSeries
    @test @inferred(a * g) isa GeneralBSeries
    # equality includes a(∅) and the trees
    @test GeneralBSeries(a) == a && GeneralBSeries(v) == v
    @test CentralBSeries(GeneralBSeries(a)) == a
    @test VectorFieldBSeries(a - b) == a - b
    @test GeneralBSeries(T, 2, coefficients(a)) != a
    @test CentralBSeries(T, zeros(Int, length(T))) != identity_series(Trees(6))
    @test isapprox(CentralBSeries(T, [0.1 * n for n in 1:length(T)]),
          CentralBSeries(T, [0.1 * n + 1e-14 for n in 1:length(T)]))
    @test CentralBSeries{Float64}(a) isa CentralBSeries{Float64} &&
          isapprox(CentralBSeries{Float64}(a), a)
    @test GeneralBSeries{Float64}(g)[empty_tree] === 5.0
    @test VectorFieldBSeries{Float64}(v) isa VectorFieldBSeries{Float64}
    @test_throws ArgumentError CentralBSeries(g)
    @test_throws ArgumentError VectorFieldBSeries(a)
    @test_throws DimensionMismatch CentralBSeries(T, [1, 2])
    @test_throws DimensionMismatch a * identity_series(Trees(4))
    @test_throws DimensionMismatch a + identity_series(Trees(4))
    # the first factor of a composition and both arguments of substitution
    # are restricted by their type
    @test_throws MethodError g * a
    @test_throws MethodError v * a
    @test_throws MethodError inv(g)
    @test_throws MethodError substitute(a, b)
    # copies are independent; slices are plain vectors
    c = copy(a); c[1] = 7
    @test a[1] != 7 && c isa CentralBSeries
    @test a[1:4] isa Vector
    # display
    s = sprint(show, MIME"text/plain"(), a)
    @test startswith(s, "CentralBSeries{Rational{BigInt}} to order 5 (17 trees)")
    @test occursin("[τ²]", s) && occursin("∅", s)
end

@testset "Ch. 3: slope, aD and implicit stages (function f)" begin
    T = Trees(6)
    D = slope(identity_series(T))
    @test D == VectorFieldBSeries(t -> order(t) == 1 ? 1 : 0, T)
    state = Ref(4)
    a = randc(T, state)
    # explicit: slope(a) = aD, (aD)([t₁ ⋯ t_m]) = a(t₁) ⋯ a(t_m)
    @test slope(a) == a * D
    @test all(slope(a)[t] == prod((a[k] for k in children(t)); init = 1)
              for t in T)
    @test slope(a, 0) == slope(a)
    # implicit: dw = wD with w = y + aii·dw
    for aii in (1 // 1, 1 // 2, big(-2) // 7)
        dw = slope(a, aii)
        @test dw == (a + aii * dw) * D
    end
    # y = id, aii = 1: implicit Euler; aii = 1/2: stage of the implicit midpoint rule
    id = identity_series(T)
    @test id + slope(id, 1) == bseries(RungeKuttaMethod(fill(1 // 1, 1, 1), [1 // 1]), T)
    @test id + slope(id, 1 // 2) ==
          bseries(RungeKuttaMethod(fill(1 // 2, 1, 1), [1 // 1]), T)
    # a diagonally implicit method, stage by stage: Fᵢ = slope(yᵢ, aᵢᵢ)
    A = [1//3 0 0; 1//2 1//4 0; -1//5 2//3 1//2]
    b = [1 // 4, 1 // 2, 1 // 4]
    F = VectorFieldBSeries[]
    for i in 1:3
        y = id + sum((A[i, j] * F[j] for j in 1:i-1); init = zero(id))
        push!(F, slope(y, A[i, i]))
    end
    @test id + sum(b[i] * F[i] for i in 1:3) == bseries(RungeKuttaMethod(A, b), T)
    @test_throws MethodError slope(D)
end

@testset "Ch. 3: λ̂ (Algorithm 6) reproduces the printed list" begin
    T = Trees(5)
    state = Ref(7)
    for _ in 1:5
        a = CentralBSeries(_ -> randq(state), T)
        @test lambda_hat(a) == book_lambda_hat(a)
        # Algorithms 5 and 6 agree
        lambda = lambda_matrix(a)
        lhat = lambda_hat(a)
        for n in 2:length(T), k in eachindex(lhat[n])
            @test lambda[n, k] == lhat[n][k]
        end
        # λ gives the same composition as the table of splittings:
        # (ab)(t) = a(t) + Σ_k λ(a, t)[k] b(t_k)
        b = GeneralBSeries(_ -> randq(state), T, randq(state))
        @test coefficients(a * b) == b[empty_tree] * coefficients(a) + lambda * coefficients(b)
    end
end

@testset "Ch. 2: splittings vs. RootedTrees.all_splittings (order 7)" begin
    T = Trees(7)
    S = splittings(T)
    for n in 1:length(T)
        ours = Dict{Any,Int}()
        for (k, F, m) in S[n]
            key = (rt(T, k), sort([rt(T, f) for f in F]))
            ours[key] = get(ours, key, 0) + m
        end
        theirs = Dict{Any,Int}()
        for (F, s) in zip(RT.all_splittings(rt(T, n))...)
            RT.order(s) == 0 && continue              # t' = ∅ is not stored
            key = (s, sort(collect(F)))
            theirs[key] = get(theirs, key, 0) + 1
        end
        @test ours == theirs
    end
end

@testset "Ch. 3: composition vs. BSeries.jl (order 6)" begin
    T = Trees(6)
    trees = [rt(t) for t in T]
    empty = RT.rootedtree(Int[])
    toBS(v) = (s = BSeries.TruncatedBSeries{eltype(trees),Rational{BigInt}}();
               s[empty] = v[empty_tree]; for n in 1:length(T); s[trees[n]] = v[n]; end; s)
    state = Ref(11)
    a = randc(T, state)
    b = randc(T, state)
    ab = compose(a, b)                         # book: first a, then b
    @test ab == a * b
    A, B = toBS(a), toBS(b)
    @test coefficients(ab) == [BSeries.compose(A, B, trees[n]) for n in 1:length(T)]
    # group properties
    E = exact_solution(T)
    @test a * inverse(a) == identity_series(T)
    @test inv(a) * a == one(a)
    c = randc(T, state)
    @test (a * b) * c == a * (b * c)
    @test a^3 == a * a * a && a^0 == one(a) && a^-2 == inv(a * a)
    # two half steps of the exact flow are one full step
    half = rescale(E, 1 // 2)
    @test half * half == E
    # 𝖡 × 𝖡* → 𝖡*: linear in the second factor, (ab)(∅) = b(∅)
    g = GeneralBSeries(_ -> big(randq(state)), T, big(3) // 2)
    v = randv(T, state)
    @test a * g isa GeneralBSeries && (a * g)[empty_tree] == 3 // 2
    @test a * v isa VectorFieldBSeries
    @test a * (2g - v) == 2 * (a * g) - a * v
    @test (a * b) * g == a * (b * g)
    @test coefficients(a * g) == [BSeries.compose(A, toBS(g), trees[n]) for n in 1:length(T)]
end

@testset "Ch. 3: substitution, modified equation (vs. BSeries.jl)" begin
    T = Trees(7)
    N = length(T)
    trees = [rt(t) for t in T]
    P = tree_partitions(T)
    # partitions = sets of cut edges: 2^(|t|-1) of them
    @test all(sum(m for (_, _, _, m) in P[n]) == 2^(order(T[n]) - 1)
              for n in 1:N)
    # pieces and skeletons agree with RootedTrees.PartitionIterator
    for n in 1:N
        ours = Dict{Any,Int}()
        for (k, F, s, m) in P[n]
            key = (sort([rt(T, f) for f in [F; k]]), rt(T, s))
            ours[key] = get(ours, key, 0) + m
        end
        theirs = Dict{Any,Int}()
        for (forest, skel) in RT.PartitionIterator(trees[n])
            key = (sort(collect(copy.(forest))), copy(skel))
            theirs[key] = get(theirs, key, 0) + 1
        end
        @test ours == theirs
    end
    empty = RT.rootedtree(Int[])
    toBS(v) = (s = BSeries.TruncatedBSeries{eltype(trees),Rational{BigInt}}();
               s[empty] = v[empty_tree]; for n in 1:N; s[trees[n]] = v[n]; end; s)
    state = Ref(5)
    b = randv(T, state)
    a = randc(T, state)
    @test substitute(b, a) isa CentralBSeries
    @test coefficients(substitute(b, a)) ==
          [BSeries.substitute(toBS(b), toBS(a), t) for t in trees]
    g = GeneralBSeries(_ -> big(randq(state)), T, -2)
    @test substitute(b, g)[empty_tree] == -2
    @test coefficients(substitute(b, g)) ==
          [BSeries.substitute(toBS(b), toBS(g), t) for t in trees]
    for rk in (euler(), runge_midpoint(), kutta4(), kutta38())
        bs = BSeries.bseries(BSeries.RungeKuttaMethod(rk.A, rk.b), 7)
        me = BSeries.modified_equation(bs)
        mi = BSeries.modifying_integrator(bs)
        @test coefficients(modified_equation(rk, T)) == [me[t] for t in trees]
        @test coefficients(modifying_integrator(rk, T)) == [mi[t] for t in trees]
    end
    # defining identities
    E = exact_solution(T)
    a = bseries(kutta4(), T)
    @test substitute(modified_equation(a), E) == a
    @test substitute(modifying_integrator(a), a) == E
    # exact flow: modified field is f itself; order p ⇒ b(t) = 0, 1 < |t| ≤ p
    @test modified_equation(E) == VectorFieldBSeries(t -> order(t) == 1 ? 1 : 0, T)
    b4 = modified_equation(a)
    @test all(b4[t] == 0 for t in T if 1 < order(t) <= 4)
    @test any(b4[t] != 0 for t in trees_of_order(T, 5))
    # explicit Euler: f − h/2 f′f + h²(1/12 f″(f,f) + 1/3 f′f′f) + …
    @test modified_equation(euler(), T)[1:4] == [1, -1 // 2, 1 // 6, 1 // 3]
end

@testset "Ch. 3, 5: Runge–Kutta methods" begin
    T = Trees(8)
    trees = [rt(t) for t in T]
    rk = kutta4()
    Phi = elementary_weights(rk, T)
    @test Phi == [RT.elementary_weight(t, rk.A, rk.b, rk.c) for t in trees]
    @test bseries(rk, T) == CentralBSeries(T, Phi)
    @test order_of_accuracy(euler(), T) == 1
    @test order_of_accuracy(runge_midpoint(), T) == 2
    @test order_of_accuracy(kutta4(), T) == 4
    @test order_of_accuracy(kutta38(), T) == 4
    @test order_of_accuracy(gauss2(), T; atol = 1e-14) == 4
    # stability function of the classical method: Σ z^k/k!, k ≤ 4
    @test stability_polynomial_coefficients(kutta4(), 4) ==
          [1, 1, 1 // 2, 1 // 6, 1 // 24]
    # RK B-series agree with BSeries.jl
    bs = BSeries.bseries(BSeries.RungeKuttaMethod(rk.A, rk.b), 6)
    for n in 1:number(last(trees_of_order(T, 6)))
        @test bs[trees[n]] == Phi[n]
    end
end

@testset "2nd edition: Butcher's seven-stage method of order 6" begin
    @test order_of_accuracy(butcher6a(), Trees(7)) == 6
end

@testset "Ch. 6: underlying one-step method of general linear methods" begin
    T = Trees(7)
    N = length(T)
    E = exact_solution(T)
    glm(m) = GeneralLinearMethod(m.A, m.U, m.B, m.V)
    methods = (BOOK_GLM_64n, BOOK_GLM_65h_i, BOOK_GLM_65h_ii, BOOK_GLM_65g,
               BOOK_GLM_64m)
    # orders claimed in the book, via Ẽ conjugate to E (Algorithm 7)
    for m in methods
        @test order_of_accuracy(glm(m), T) == m.order
    end
    @test is_partitioned_form(glm(BOOK_GLM_64n))
    @test !is_partitioned_form(glm(BOOK_GLM_65g))
    # an RK method as a GLM (r = 1): Ẽ is the RK method itself
    for rk in (kutta4(), kutta38(), butcher6a())
        res = underlying_one_step(GeneralLinearMethod(rk), T)
        @test res.E == bseries(rk, T)
        @test res.S == [identity_series(T)]
    end
    # a perturbed coefficient destroys the order
    m = BOOK_GLM_65h_ii
    B = copy(m.B); B[1, 4] += 1 // 1000; B[1, 1] -= 1 // 1000
    @test order_of_accuracy(GeneralLinearMethod(m.A, m.U, B, m.V), T) < 5
    tall = [T[1]]                               # τ, [τ], [[τ]], …
    for _ in 2:7
        push!(tall, bplus(T, [tall[end]]))
    end
    for m in methods
        g = glm(m)
        res = underlying_one_step(g, T)
        (; A, U, V) = g
        Bm = g.B
        r, s = size(V, 1), size(A, 1)
        # M ∘ S = S ∘ Ẽ in the original input basis, all components
        etaD = Vector{Vector{eltype(V)}}(undef, N)
        eta = similar(etaD)
        for n in 1:N
            etaD[n] = n == 1 ? ones(eltype(V), s) :
                    etaD[T.table.L[n]] .* eta[T.table.R[n]]
            eta[n] = A * etaD[n] + U * [res.S[i][n] for i in 1:r]
        end
        out = [Bm * etaD[n] + V * [res.S[i][n] for i in 1:r] for n in 1:N]
        for i in 1:r
            @test res.S[i][empty_tree] == res.u[i]
            @test coefficients(res.E * res.S[i]) == [out[n][i] for n in 1:N]
        end
        @test V * res.u == res.u && U * res.u == ones(s)
        # linear problems: Ẽ on the tall trees is the principal
        # eigenvalue w(z) = 1 + Σ w_k z^k of M(z) = V + Σ z^k B A^(k-1) U,
        # computed by perturbation theory without trees
        Id = [i == j ? 1 // 1 : 0 // 1 for i in 1:r, j in 1:r]
        u = res.u
        wl = BSeriesBook.solve_exact([(V - Id)'; u'], [zeros(eltype(V), r); 1])
        Mk = [Bm * A^(k - 1) * U for k in 1:7]
        v = [u]
        w = [big(1) // 1]
        for k in 1:7
            wk = sum(wl' * Mk[j] * v[k-j+1] for j in 1:k)
            rhs = sum((w[j+1] * v[k-j+1] for j in 1:k-1); init = zeros(eltype(V), r)) -
                  sum(Mk[j] * v[k-j+1] for j in 1:k) + wk * u
            push!(v, BSeriesBook.solve_exact([V - Id; wl'], [rhs; 0]))
            push!(w, wk)
        end
        @test [res.E[tall[k]] for k in 1:7] == w[2:8]
    end
end

@testset "Ch. 6: conjugacy (Algorithm 7, (6.6 e))" begin
    T = Trees(5)
    state = Ref(3)
    for _ in 1:5
        b = randc(T, state)
        b[1] = 1
        x = randc(T, state)
        a = inv(x) * b * x
        # (6.6 e): the order-4 conjugacy invariants of a and b agree
        @test book_conjugate4(a) == book_conjugate4(b)
    end
    # One step of the order-by-order conjugation: a agrees with b for
    # |t| < p, and a = x⁻¹ b x with x supported on order p − 1.
    for p in 2:5, _ in 1:3
        b = randc(T, state)
        b[1] = 1
        x = identity_series(T)
        for t in trees_of_order(T, p - 1)
            x[t] = randq(state)
        end
        a = inv(x) * b * x
        @test all(a[t] == b[t] for t in T if order(t) < p)
        res = conjugacy_step(a, b, p)
        @test res !== nothing
        @test all(res.conjugated[t] == b[t] for t in T if order(t) <= p)
        @test res.conjugated == res.x * a * inv(res.x)
    end
    # Order 4: a and b agree below order 4 but violate (6.6 e) — not conjugate
    for _ in 1:5
        b = randc(T, state)
        b[1] = 1
        a = copy(b)
        a[8] += 1                     # a₈ ≠ b₈
        @test conjugacy_step(a, b, 4) === nothing
    end
end

end
