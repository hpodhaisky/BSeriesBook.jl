using Test
using BSeriesBook
import BSeries
import RootedTrees as RT

include("book_values.jl")

# Map from the standard numbering to RootedTrees.jl trees (level sequences)
rt(T, n) = RT.rootedtree(level_sequence(T, n))

randq(rng_state) = (rng_state[] = (1103515245 * rng_state[] + 12345) % 2^31;
                    (rng_state[] % 19 - 9) // (rng_state[] % 7 + 1))

@testset "BSeriesBook" begin

@testset "Ch. 2: Algorithm 3, generate trees (book values to order 5)" begin
    T = generate_trees(5)
    @test length(T) == 17
    @test T.L[2:17] == BOOK_L
    @test T.R[2:17] == BOOK_R
    for (ℓ, row) in BOOK_prod, (r, n) in enumerate(row)
        @test T.prod[ℓ, r] == n
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
        t = rt(T, n)
        push!(seen, t)
        @test RT.order(t) == order(T, n)
        @test RT.symmetry(t) == σ(T, n)
        @test RT.density(t) == γ(T, n)
        @test RT.α(t) == α(T, n)
        @test RT.β(t) == β(T, n)
        @test tree_number(T, children(T, n)) == n
    end
    @test length(seen) == length(T)            # all trees distinct
    # β-product: t_ℓ ∘ t_r is the Butcher product `∘` in RootedTrees.jl
    for ℓ in 1:length(T), r in 1:length(T)
        order(T, ℓ) + order(T, r) ≤ 8 || continue
        @test rt(T, beta_product(T, ℓ, r)) ==
              rt(T, ℓ) ∘ rt(T, r)
    end
end

@testset "Ch. 3: λ̂ (Algorithm 6) reproduces the printed list" begin
    T = Trees(5)
    state = Ref(7)
    for _ in 1:5
        a = [randq(state) for _ in 1:length(T)]
        @test lambda_hat(T, a, 5) == book_lambda_hat(a)
        # Algorithms 5 and 6 agree
        λ = lambda_matrix(T, a, 5)
        λ̂ = lambda_hat(T, a, 5)
        for n in 2:length(T), k in eachindex(λ̂[n])
            @test λ[n, k] == λ̂[n][k]
        end
        # λ gives the same composition as the table of splittings:
        # (ab)(t) = a(t) + Σ_k λ(a, t)[k] b(t_k)
        b = [randq(state) for _ in 1:length(T)]
        @test a .+ λ * b == compose(a, b, T)
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
    trees = [rt(T, n) for n in 1:length(T)]
    empty = RT.rootedtree(Int[])
    toBS(v) = (s = BSeries.TruncatedBSeries{eltype(trees),Rational{BigInt}}();
               s[empty] = 1; for n in 1:length(T); s[trees[n]] = v[n]; end; s)
    state = Ref(11)
    a = [big(randq(state)) for _ in 1:length(T)]
    b = [big(randq(state)) for _ in 1:length(T)]
    ab = compose(a, b, T)                      # book: first a, then b
    A, B = toBS(a), toBS(b)
    @test ab == [BSeries.compose(A, B, trees[n]) for n in 1:length(T)]
    # group properties
    E = exact_solution(T)
    @test compose(a, inverse(a, T), T) == identity_series(T)
    @test compose(inverse(a, T), a, T) == identity_series(T)
    c = [big(randq(state)) for _ in 1:length(T)]
    @test compose(compose(a, b, T), c, T) == compose(a, compose(b, c, T), T)
    # two half steps of the exact flow are one full step
    half = [E[n] / 2^order(T, n) for n in 1:length(T)]
    @test compose(half, half, T) == E
end

@testset "Ch. 3: substitution, modified equation (vs. BSeries.jl)" begin
    T = Trees(7)
    N = length(T)
    trees = [rt(T, n) for n in 1:N]
    P = tree_partitions(T)
    # partitions = sets of cut edges: 2^(|t|-1) of them
    @test all(sum(m for (_, _, _, m) in P[n]) == 2^(order(T, n) - 1)
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
    toBS(v, v0) = (s = BSeries.TruncatedBSeries{eltype(trees),Rational{BigInt}}();
                   s[empty] = v0; for n in 1:N; s[trees[n]] = v[n]; end; s)
    state = Ref(5)
    b = [big(randq(state)) for _ in 1:N]
    a = [big(randq(state)) for _ in 1:N]
    @test substitute(b, a, T) ==
          [BSeries.substitute(toBS(b, 0), toBS(a, 1), t) for t in trees]
    for rk in (euler(), runge_midpoint(), kutta4(), kutta38())
        bs = BSeries.bseries(BSeries.RungeKuttaMethod(rk.A, rk.b), 7)
        me = BSeries.modified_equation(bs)
        mi = BSeries.modifying_integrator(bs)
        @test modified_equation(rk, T) == [me[t] for t in trees]
        @test modifying_integrator(rk, T) == [mi[t] for t in trees]
    end
    # defining identities
    E = exact_solution(T)
    a = bseries(kutta4(), T)
    @test substitute(modified_equation(a, T), E, T) == a
    @test substitute(modifying_integrator(a, T), a, T) == E
    # exact flow: modified field is f itself; order p ⇒ b(t) = 0, 1 < |t| ≤ p
    @test modified_equation(E, T) == [n == 1 ? 1 : 0 for n in 1:N]
    b4 = modified_equation(a, T)
    @test all(b4[n] == 0 for n in 2:N if order(T, n) ≤ 4)
    @test any(b4[n] != 0 for n in trees_of_order(T, 5))
    # explicit Euler: f − h/2 f′f + h²(1/12 f″(f,f) + 1/3 f′f′f) + …
    @test modified_equation(euler(), T)[1:4] == [1, -1 // 2, 1 // 6, 1 // 3]
end

@testset "Ch. 3, 5: Runge–Kutta methods" begin
    T = Trees(8)
    trees = [rt(T, n) for n in 1:length(T)]
    rk = kutta4()
    Φ = elementary_weights(rk, T)
    @test Φ == [RT.elementary_weight(t, rk.A, rk.b, rk.c) for t in trees]
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
    for n in 1:last(trees_of_order(T, 6))
        @test bs[trees[n]] == Φ[n]
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
        @test underlying_one_step(GeneralLinearMethod(rk), T).E ==
              bseries(rk, T)
    end
    # a perturbed coefficient destroys the order
    m = BOOK_GLM_65h_ii
    B = copy(m.B); B[1, 4] += 1 // 1000; B[1, 1] -= 1 // 1000
    @test order_of_accuracy(GeneralLinearMethod(m.A, m.U, B, m.V), T) < 5
    tall = [1]                                  # τ, [τ], [[τ]], …
    for _ in 2:7
        push!(tall, tree_number(T, [tall[end]]))
    end
    for m in methods
        g = glm(m)
        res = underlying_one_step(g, T)
        (; A, U, V) = g
        Bm = g.B
        r, s = size(V, 1), size(A, 1)
        # M ∘ S = S ∘ Ẽ in the original input basis, all components
        ηD = Vector{Vector{eltype(V)}}(undef, N)
        η = similar(ηD)
        for n in 1:N
            ηD[n] = n == 1 ? ones(eltype(V), s) :
                    ηD[T.table.L[n]] .* η[T.table.R[n]]
            η[n] = A * ηD[n] + U * res.ζ[n]
        end
        out = [Bm * ηD[n] + V * res.ζ[n] for n in 1:N]
        for i in 1:r
            ζi = [res.ζ[n][i] for n in 1:N]
            @test compose(res.E, ζi, T; b0 = res.u[i]) ==
                  [out[n][i] for n in 1:N]
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
        b = [big(randq(state)) for _ in 1:length(T)]
        b[1] = 1
        x = [big(randq(state)) for _ in 1:length(T)]
        a = compose(compose(inverse(x, T), b, T), x, T)   # a = x⁻¹ b x
        # (6.6 e): the order-4 conjugacy invariants of a and b agree
        @test book_conjugate4(a) == book_conjugate4(b)
    end
    # One step of the order-by-order conjugation: a agrees with b for
    # |t| < p, and a = x⁻¹ b x with x supported on order p − 1.
    for p in 2:5, _ in 1:3
        b = [big(randq(state)) for _ in 1:length(T)]
        b[1] = 1
        x = zero(b)
        for t in trees_of_order(T, p - 1)
            x[t] = randq(state)
        end
        a = compose(compose(inverse(x, T), b, T), x, T)
        @test all(a[t] == b[t] for t in 1:length(T) if order(T, t) < p)
        ok, a′, _ = conjugacy_step(a, b, T, p)
        @test ok
        @test all(a′[t] == b[t] for t in 1:length(T) if order(T, t) ≤ p)
    end
    # Order 4: a and b agree below order 4 but violate (6.6 e) — not conjugate
    for _ in 1:5
        b = [big(randq(state)) for _ in 1:length(T)]
        b[1] = 1
        a = copy(b)
        a[8] += 1                     # a₈ ≠ b₈
        @test conjugacy_step(a, b, T, 4) == (false,)
    end
end

end
