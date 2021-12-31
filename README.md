# AdaptiveWindow.jl
Adaptive window implementation (ADWIN2: Bifet &amp; Gavalda) 


# Example

```julia
using AdaptiveWindow

m = AdaptiveMean(0.001) 

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
The package uses the interface from OnlineStatsBase. The 'nobs' measurement returns the 
number of observations underlying the current mean value. Values might have dropped out of the
adaptive window.

The fit! method returns the statistic itself. In some cases, when you want to detect windows being dropped. To 

the update! routine does the work and returns True if values are dropped.

In the following example, the adaptive window should detect the shift in distribution about 50 observations after the 500th one. If you want the mean to be more sensitive increase the probability of detecting shift (delta). Note that this will lead more often to false positives

```julia
using AdaptiveWindow

m = AdaptiveMean(0.001, ad -> println("Shift detected!"))

for i in 1:1_000
    r = randn()
    if i > 500
        r += 1
    end
    fit!(m, r)
end

```




