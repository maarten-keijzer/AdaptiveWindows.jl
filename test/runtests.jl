# Adaptive

using AdaptiveWindow
using Test

ad = AdWin(0.001)

# This should not trigger a truncated window
fit!(ad, randn(10_000))
@test stats(ad).n == 10_000

# This should trigger a truncated window
fit!(ad, 1 .+ randn(10_000))
@test 9_900 < stats(ad).n < 20_000


