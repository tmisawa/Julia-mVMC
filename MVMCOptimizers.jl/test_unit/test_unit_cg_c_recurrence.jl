using Test, LinearAlgebra
using MVMCOptimizers

@testset "C 1.4.0 CG stops a tiny SPD system at the absolute 1e-20 guard" begin
    fields = split(read(joinpath(@__DIR__, "fixtures", "c_cg_small_spd.txt"), String))
    status = parse(Int, fields[1])
    solution, residual, direction = parse.(Float64, fields[2:4])
    # Binary-exact one-dimensional SPD: O=2^-54, S=2^-108, g=2^-54.
    # delta=2^-108 is above the threshold 1e-60, but d*S*d=2^-216 < 1e-20, so the
    # C 1.4.0 Main returns -(0+1) before updating x (1.3.0 returned x=2^54).
    @test status == -1
    for complex in (false, true)
        for diag_scale in (false, true)
            ws = MVMCOptimizers.CGWorkspace(1, 1, complex)
            ws.stcOs_real[1,1] = ldexp(1.0, -54)
            ws.g[1] = ldexp(1.0, -54)
            ws.sdiag[1] = 1.0  # only read by DiagScale; keeps d = r there
            main! = diag_scale ? MVMCOptimizers.stochastic_opt_precon_cg_diag_scale_main! :
                                 MVMCOptimizers.stochastic_opt_cg_main!
            result = main!(ws, 1, 1, 1.0, 0.0, 1e-30, 1, complex)
            @test result == status
            @test ws.x[1] == solution  # untouched initial zero
            @test abs(ws.r[1] - residual) <= eps(Float64) * ldexp(1.0, -54)
            diag_scale || @test abs(ws.d[1] - direction) <= eps(Float64) * ldexp(1.0, -54)
        end
    end
end

@testset "C 1.4.0 sampled operator and CG residuals across limits 1 through 41" begin
    decode(line) = [reinterpret(Float64, parse(UInt64, word; base=16)) for word in split(line)]
    for name in ("real", "complex", "sampled_complex",
                 "real_diag_scale", "complex_diag_scale", "sampled_complex_diag_scale")
        diag_scale = endswith(name, "_diag_scale")
        main! = diag_scale ? MVMCOptimizers.stochastic_opt_precon_cg_diag_scale_main! :
                             MVMCOptimizers.stochastic_opt_cg_main!
        lines = filter(line -> !startswith(line, "#"), readlines(joinpath(@__DIR__, "fixtures", "c_cg_refresh", name * ".txt")))
        @test length(lines) == 7 + 3*41
        n, samples, complex_int = parse.(Int, split(lines[1]))
        @test n > 0
        @test samples > 0
        @test complex_int in (0, 1)
        complex = complex_int != 0
        mean, diagonal, real_samples, imag_samples, gradient, expected_product = decode.(lines[2:7])
        @test length(mean) == n
        @test length(diagonal) == n
        @test length(gradient) == n
        @test length(expected_product) == n
        @test length(real_samples) == n*samples
        @test length(imag_samples) == (complex ? n*samples : 0)
        function workspace()
            ws = MVMCOptimizers.CGWorkspace(n, samples, complex)
            ws.stcO .= mean; ws.sdiag .= diagonal; ws.g .= gradient
            ws.stcOs_real .= reshape(real_samples, n, samples)
            complex && (ws.stcOs_imag .= reshape(imag_samples, n, samples))
            return ws
        end
        ws = workspace()
        product = zeros(n)
        MVMCOptimizers.operate_by_s!(product, gradient, ws, n, samples, 1.0/samples, 1e-5, complex)
        operator_budget = 8 * (n + samples) * eps(Float64)
        @test all(abs(a-b) <= operator_budget + operator_budget*abs(b) for (a,b) in zip(product, expected_product))
        # Independent explicit covariance; does not invoke the sampled operator.
        covariance = zeros(n,n)
        for row in 1:n, col in 1:n
            gram = 0.0
            for sample in 1:samples
                r = (sample-1)*n + row; c = (sample-1)*n + col
                gram += real_samples[r]*real_samples[c] + (complex ? imag_samples[r]*imag_samples[c] : 0.0)
            end
            covariance[row,col] = gram/samples - mean[row]*mean[col] + (row==col ? 1e-5*diagonal[row] : 0.0)
        end
        for limit in 1:41
            base = 8 + 3*(limit-1)
            words = split(lines[base]; limit=3)
            @test parse(Int,words[1]) == limit
            # C 1.4.0 status: tolerance 0 never converges, so every limit is
            # -(limit+1), or -(iter+1) at an earlier breakdown guard.
            expected_status = parse(Int,words[2])
            @test -(limit+1) <= expected_status < 0
            expected = (decode(words[3]), decode(lines[base+1]), decode(lines[base+2]))
            @test all(length(vector) == n for vector in expected)
            ws = workspace()
            status = main!(ws,n,samples,1.0/samples,1e-5,0.0,limit,complex)
            @test status == expected_status
            # These archived problems have condition numbers of 1e7--1e9.
            # Later truncated iterates are sensitive to BLAS reduction order;
            # iteration * eps is not a cross-platform forward-error bound.
            # Keep the C first-step comparison and validate BOTH trajectories
            # against the independent covariance at every later limit.
            budget = 8*(n+samples)*eps(Float64)
            # DiagScale: r_new = r - alpha*q cancels (rounding ~ eps*(|g|+|r_new|)),
            # and d_new = M^-1 r_new + beta*d divides it by M = (1+shift)*sdiag,
            # which is as small as 1e-12. Bound those components by the operand
            # scales they are computed from instead of by their own magnitude.
            precon = [abs(MVMCOptimizers._cg_diag_scale(diagonal[i], 1 + 1e-5)) for i in 1:n]
            scale_r = diag_scale ? abs.(gradient) .+ abs.(expected[2]) : zeros(n)
            scale_d = diag_scale ? (2 .* abs.(gradient) .+ abs.(expected[2])) ./ precon : zeros(n)
            for (k, (actual,target)) in enumerate(zip((ws.x,ws.r,ws.d),expected))
                @test length(actual) == length(target)
                @test all(isfinite, actual)
                @test all(isfinite, target)
                if limit == 1
                    extra = k == 2 ? scale_r : k == 3 ? scale_d : zeros(n)
                    @test all(abs(actual[i]-target[i]) <= budget*(1+abs(target[i])+extra[i])
                              for i in 1:n)
                end
            end
            for (solution, residual) in ((ws.x, ws.r), (expected[1], expected[2]))
                for row in 1:n
                    ax = 0.0; scale = abs(gradient[row])
                    for col in 1:n
                        term = covariance[row,col]*solution[col]
                        ax += term; scale += abs(term)
                    end
                    residual_budget = 16*(n+samples)*limit*eps(Float64)*scale
                    @test abs(residual[row]-(gradient[row]-ax)) <= residual_budget
                end
            end
            if limit % 20 == 0 && status == -(limit+1)
                # A refresh must replace the accumulated residual with g-S*x.
                # Identical inputs on the same backend are deterministic; this
                # is not a bitwise comparison against a different C provider.
                refreshed = zeros(n)
                MVMCOptimizers.operate_by_s!(refreshed,ws.x,ws,n,samples,1.0/samples,1e-5,complex)
                @test ws.r == gradient - refreshed
            end
        end
        @test length(lines) == 7 + 3*41
    end
end

@testset "sampled CG converges to an independent analytic SPD solution" begin
    # Paired samples m +/- a_i*e_i have mean m and covariance diag(a_i^2/n).
    # Imaginary samples +/- b_i*e_i add diag(b_i^2/n) with zero imaginary mean.
    # Including the diagonal shift gives a known spectrum and solution without
    # reproducing the CG recurrence or calling the production matrix operator.
    n = 64
    samples = 2n
    shift = 0.125
    mean = [(-1.0)^i / 64 for i in 1:n]
    a = collect(1:n) ./ 8
    # C 1.4.0 stops at the absolute guards |d.q| < 1e-20 and old |delta| < 1e-20.
    # Scaling the exact solution by 2^27 keeps the relative stop threshold
    # (rtol*|g|)^2 ~ 1e-15 and every d.q above them, so the relative convergence
    # gate is reached before any absolute guard.
    exact = [ldexp((-1.0)^i / (i+1), 27) for i in 1:n]
    for complex in (false, true), diag_scale in (false, true)
        b = complex ? a ./ 2 : zeros(n)
        diagonal = (a.^2 + b.^2) ./ n
        eigenvalues = (1 + shift) .* diagonal
        ws = MVMCOptimizers.CGWorkspace(n, samples, complex)
        ws.stcO .= mean
        ws.sdiag .= diagonal
        ws.g .= eigenvalues .* exact
        ws.stcOs_real .= mean
        for i in 1:n
            ws.stcOs_real[i,2i-1] += a[i]
            ws.stcOs_real[i,2i] -= a[i]
            if complex
                ws.stcOs_imag[i,2i-1] = b[i]
                ws.stcOs_imag[i,2i] = -b[i]
            end
        end
        rtol = 1e-14
        tolerance = rtol * norm(ws.g) / n
        main! = diag_scale ? MVMCOptimizers.stochastic_opt_precon_cg_diag_scale_main! :
                             MVMCOptimizers.stochastic_opt_cg_main!
        iterations = main!(ws,n,samples,1.0/samples,shift,tolerance,10n,complex)
        if diag_scale
            # The Jacobi preconditioner (1+shift)*sdiag equals this diagonal S.
            @test 0 <= iterations <= 2
        else
            @test 40 < iterations < 10n # Both residual refresh boundaries are crossed.
        end
        @test norm(ws.r) < rtol * norm(ws.g)
        condition = maximum(eigenvalues) / minimum(eigenvalues)
        @test condition == 4096
        # Dimension- and condition-scaled roundoff bounds, independent of the
        # C fixtures. The backward check also detects a wrong mean/imag/shift.
        budget = 32n * eps(Float64)
        @test norm(ws.x - exact) <= budget * condition * norm(exact)
        @test norm(ws.g - eigenvalues .* ws.x) <=
            budget * (norm(ws.g) + maximum(eigenvalues)*norm(ws.x))
    end
end
