
using AdaptiveWindow
using Test

@testset "Mean Computation" begin
    m = AdaptiveMean(1e-9)

    r = randn(1000)

    fit!(m, r)

    m1 = sum(r) / length(r)
    m2 = value(m)

    @test m1 ≈ m2

end

@testset "Distribution Shift" begin
    ad = AdaptiveMean(0.001)

    # This should not trigger a truncated window
    fit!(ad, randn(10_000))
    @test stats(ad).n == 10_000

    # Changing the distribution should trigger a truncated window
    fit!(ad, 1 .+ randn(10_000))
    @test 9_900 < stats(ad).n < 20_000

    # check truncation of shifting using the callback function
    shifted = false

    m = AdaptiveMean(0.001, onshiftdetected = ad -> shifted = true)

    for i in 1:1_000
        r = randn()
        if i > 500
            r += 1
        end
        fit!(m, r)
    end

    @test shifted
end

@testset "Memory Usage" begin
    m = AdaptiveMean()
    n = 100_000

    # Run without dropping observations -- for speed
    fit!(withoutdropping(m), randn(n))

    @test length(m.window) < AdaptiveWindow.M * log2(n)
end
