# AdaptiveWindow.jl
Adaptive window implementation (ADWIN2: Bifet &amp; Gavalda) 

The adaptive windowing algorithm works by maintaining a set of buckets (the window) with Variance() objects in each of these. Each group of buckets has a maximum capacity for the number of observations it holds. A constant M is used to define the size of the group. In the package, M is hardcoded to be equal to 5, but M=3 will be used below.

The content of the window is exponentially increasing, with the first group of buckets containing at most 1 observation, the next group 2, then 4, etc.:

1 1 1 2 2 2 4 4 4 8 8 8 16 16 16 32 32 32 ...

The algorithm perform a statistical test using Hoeffding bounds on every splitpoint between these buckets, using the accumulated means and variances. If the front of the window is considered to have a statistically different mean than the back of the window, the observations at the back are removed, and a 'drift' is detected. 


# Example

```julia
using AdaptiveWindow

m = AdaptiveMean(δ = 0.001) 

fit!(m, randn(1_000))

println(nobs(m)) # we should see that the stats are computed over 1_000 data points

prinln("The mean: $(value(m)) should be close to 0.0")

# change the distribution
fit!(m, randn(1_000) .+ 1) 

# now the older observations should be dropped
println(nobs(m)) 

# mean should be close to 1.0 (not 0.5 if we did a mean over all 2_000 points)
prinln("The mean: $(value(m)) should be close to 1.0")

```
The package uses the interface from OnlineStatsBase. The 'nobs' measurement returns the number of observations underlying the current mean value. Values might have dropped out of the adaptive window, so nobs is not necessarily equal to the number of times fit! is called.

In the following example, the adaptive window should detect the shift in distribution about 50 observations after the 500th one. If you want the mean to be more sensitive increase the probability of detecting shift (delta). Note that this will lead more often to false positives.

```julia
using AdaptiveWindow

m = AdaptiveMean(onshiftdetected = ad -> println("Shift detected!"))

for i in 1:1_000
    r = randn()
    if i > 500
        r += 1
    end
    fit!(m, r)
end

```




