
using AdaptiveWindows
using Test
using Distributions

@testset verbose=true "Adaptive Mean" begin

    @testset "Mean Computation " begin
        m = AdaptiveMean(δ = 1e-9)

        r = randn(1000)

        fit!(m, r)

        m1 = sum(r) / length(r)
        m2 = value(m)

        @test m1 ≈ m2
        ad = AdaptiveMean()

        # This should not trigger a truncated window
        fit!(ad, randn(10_000))
        @test stats(ad).n == 10_000

        # Changing the distribution should trigger a truncated window
        fit!(ad, 1 .+ randn(10_000))
        @test 9_900 < stats(ad).n < 20_000

        # check truncation of shifting using the callback function
        shifted = false

        m = AdaptiveMean(onshiftdetected = (ad, idx) -> shifted = true)

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

        mn = AdaptiveMean()
        n = 1<<12

        # withoutdropping for speed
        fit!(withoutdropping(mn), ones(n))
        m = AdaptiveWindows.M 
        expected = m * ceil(log2(n) - log2(m))
        @test length(mn.window) == expected
        @test nobs(mn) == n 
        @test consistent(mn)
        
        # Maximum amount of memory
        mn = withmaxlength(AdaptiveMean(), 3)
        fit!(mn, rand(10000))
        @test length(mn.ad.window) == AdaptiveWindows.M * 3
        @test consistent(mn.ad)
    end
    
end

@testset verbose=true "Synced AdWin" begin
    m = AdWinGroup(3)
    fit!(m, [1, 0, 0])
    @test nobs(m) == 1
    fit!(m, [0, 1, 0])
    @test nobs(m) == 2
    fit!(m, [0, 0, 1])
    @test nobs(m) == 3
    fit!(m, [1, 0, 0])
    @test nobs(m) == 4

    @test mean.(m.adwins) == [0.5, 0.25, 0.25]

    m = AdWinGroup(3)
    #start with uniform distribution
    for i in 1:999
        vals = zeros(3)
        vals[mod1(i, 3)] = 1
        fit!(m, vals)
    end
    @test mean.(m.adwins) ≈ [1/3, 1/3, 1/3]
    @test nobs(m) == 999
    # now shift to a non-uniform distribution
    for i in 1:10000
        vals = zeros(3)
        vals[mod1(i % 4, 3)] = 1
        fit!(m, vals)
        @test mean.(m.adwins) |> sum ≈ 1
    end
    # check if it shifted
    @test nobs(m) <= 10000
    for adwin in m.adwins
        @test nobs(adwin) == nobs(m)
    end
    @test mean.(m.adwins) |> sum ≈ 1

    m = AdWinGroup(2)

    for i in 1:20000
        prob = [sin(i/10000 * 2π)^2, cos(i/10000 * 2π)^2]
        prob = prob ./ sum(prob)
        class = rand(Categorical(prob))
        res = zeros(2)
        res[class] = 1.0
        fit!(m, res)

        @test mean.(m.adwins) |> sum ≈ 1

    end

end

