
using AdaptiveWindows
using Test

@testset verbose=true "Adaptive Mean" begin

    @testset "Mean Computation" begin
        m = AdaptiveMean(δ = 1e-9)

        r = randn(1000)

        fit!(m, r)

        m1 = sum(r) / length(r)
        m2 = value(m)

        @test m1 ≈ m2

    end

    @testset "Distribution Shift" begin
        ad = AdaptiveMean()

        # This should not trigger a truncated window
        fit!(ad, randn(10_000))
        @test stats(ad).n == 10_000

        # Changing the distribution should trigger a truncated window
        fit!(ad, 1 .+ randn(10_000))
        @test 9_900 < stats(ad).n < 20_000

        # check truncation of shifting using the callback function
        shifted = false

        m = AdaptiveMean(onshiftdetected = ad -> shifted = true)

        for i in 1:1_000
            r = randn()
            if i > 500
                r += 1
            end
            fit!(m, r)
        end

        @test shifted
    end

    function consistent(ad)
        total = sum(nobs(v) for v in ad.window)
        total == nobs(ad.stats)
    end

    @testset "Memory Management" begin

        m = AdaptiveMean()
        fit!(m, 1)
        @test nobs(m.window[1]) == 0
        @test nobs(m.window[2]) == 1
        @test nobs(m.window[3]) == 0
        fit!(m, 1)
        @test nobs(m.window[1]) == 0
        @test nobs(m.window[2]) == 1
        @test nobs(m.window[3]) == 1
        fit!(m, 1)
        fit!(m, 1)
        fit!(m, 1)
        fit!(m, 1)
        fit!(m, 1)
        @test consistent(m)
    
        m = AdaptiveMean()
        n = AdaptiveWindows.M * ( 1 + 2 + 4)
        fit!(m, ones(n))

        @test length(m.window) <= AdaptiveWindows.M * log2(n)
        @test nobs(m) == n 
        @test consistent(m)

        m = AdaptiveMean()
        n = 10_000

        # withoutdropping for speed
        fit!(withoutdropping(m), ones(n))
        @test length(m.window) <= AdaptiveWindows.M * log2(n)
        @test nobs(m) == n 
        @test consistent(m)
    end
end
