module AdaptiveWindows

export AdaptiveMean, fit!, value, mean, nobs, stats, withoutdropping, withmaxlength

import StatsBase: nobs, fit!, merge!
import OnlineStatsBase: value, OnlineStat, Variance, Mean, _fit!

#=
 Adaptive Windowing version 2 (AdaptiveMean2)

    Used to track the mean value of a stream of data with a possibly changing population.
    If the population is determined to be changed, older observations will be dropped

  Bifet and Gavalda. Learning from Time-Changing Data with Adaptive Windowing
=#

    # Bifet and Gavalda: We use, somewhat arbitrarily, M = 5 for all experiments.
    const M = 5

    abstract type AdaptiveBase <: OnlineStat{Number} end

    """
        Tracks an "adaptive" mean by keeping a vector of variance objects of exponentially
        increating number of observations, performing a test based on Hoeffding bounds that 
        assesses if older events are drawn from a different distribution than new events,
        i.e., did the distribution shift?
    """
    mutable struct AdaptiveMean <: AdaptiveBase
        
        δ ::Float64
        window::Vector{Variance}
        stats ::Variance
        
        onshiftdetect

        AdaptiveMean(;δ = 0.001, onshiftdetected = identity) = new(δ, [Variance() for _ in 1:M], Variance(), onshiftdetected)    
    end

    function _fit!(ad::AdaptiveMean, value)
        fit!(ad.window[1], value)
        fit!(ad.stats, value)
    
        compress!(ad)
        if dropifdrifting!(ad)
            ad.onshiftdetect(ad)
        end
        ad
    end

    nobs(ad::AdaptiveMean) = ad.stats.n
    stats(ad::AdaptiveMean) = ad.stats

    function mean(ad::AdaptiveMean)
        ad.stats.μ
    end

    function value(ad::AdaptiveMean)
        mean(ad)
    end
                
    function compress!(ad::AdaptiveMean)
        makespace!(ad, 1, 1.0);
    end
    
    function print_trace(m)
        for v in m.window
            print(nobs(v), " ")
        end
        println("stats: ", nobs(m.stats))
    end    

    function makespace!(ad::AdaptiveMean, start::Int, max::Float64)
    #=
        The window is a gappy list of data points, this avoids allocations and reallocations 
        when the window resizes
    
        Assume M = 3, the lists below show the counts of observations
    
            Entry: [1, 0, 0],               Exit: [0, 1, 0]
            Entry: [1, 1, 0],               Exit: [0, 1, 1]
            Entry: [1, 1, 1],               Exit: [0, 1, 1, 1, 0, 0]
            Entry [1, 1, 1, 1, 0, 0],       Exit: [0, 1, 1, 0, 2, 0]
            Entry [1, 1, 1, 0, 2, 0],       Exit: [0, 1, 1, 1, 2, 0]
            Entry [1, 1, 1, 1, 2, 0],       Exit: [0, 1, 1, 0, 2, 2]
            Entry [1, 1, 1, 0, 2, 2],       Exit: [0, 1, 1, 1, 2, 2]
            Entry [1, 1, 1, 1, 2, 2],       Exit: [0, 1, 1, 0, 2, 2, 2, 0, 0]
    
    =#
    
        if nobs(ad.window[start]) < max
            return
        end
    
        # Move-to-front: [a, b, c, d] -> [d, a, b, c]
        last_entry_in_range = start - 1 + M;
        laststats = ad.window[last_entry_in_range];
        for j in last_entry_in_range:-1:start+1
            ad.window[j] = ad.window[j-1];
        end
        ad.window[start] = laststats

        if nobs(ad.window[start]) != 0
            next = start + M
            if length(ad.window) < next
                for i in 1:M
                    push!(ad.window, Variance())
                end
            end
    
            merge!(ad.window[next], laststats)
            ad.window[start] = Variance()
            makespace!(ad, next, max*2)
        end
    end

    function tomean(v::Variance)
        mean = Mean()
        mean.n = v.n 
        mean.μ = v.μ
        mean 
    end 

    function _remove!(m::Mean, v::Variance)
        residualsum = m.μ * m.n - v.μ * v.n
        m.n -= v.n 
        m.μ = residualsum / m.n 
        m
    end
       
    function _merge!(m::Mean, v::Variance)
        merge!(m, tomean(v))
    end

    function dropifdrifting!(ad::AdaptiveMean)
    
        statsToRight = tomean(ad.stats)
        statsToLeft = Mean()
    
        deltaPrime = ad.δ / log(nobs(ad.stats))
        logDeltaPrime = log(2/deltaPrime)
        variance = value(ad.stats)
    
        for i in 2:length(ad.window)
            if nobs(ad.window[i]) == 0
                continue
            end
    
            _remove!(statsToRight, ad.window[i])
            _merge!(statsToLeft, ad.window[i])
    
            if statsToRight.n < 1e-9
                break # only zeros from this point on
            end
    
            mInv = 1.0/ statsToRight.n + 1.0/ statsToLeft.n;
    
            epsCut =sqrt(2 * mInv * variance * logDeltaPrime) + 2.0/3.0 * mInv * logDeltaPrime;
    
            if abs(statsToRight.μ - statsToLeft.μ) > epsCut
                # Drift detected, clear all from here on
                for j = i+1:length(ad.window)
                    ad.window[j] = Variance()
                end
                ad.stats = Variance()
                for stat in 1:i
                    merge!(ad.stats, ad.window[stat])
                end
                return true
            end
        end
    
        return false
    end
    
    struct NoDropWrapper <: AdaptiveBase
        ad::AdaptiveBase
    end

    """
        Creates a variant of an AdaptiveMean that updates (fit!) an AdaptiveMean and prevents
        statistical tests being performed while updating the data to improve performance. 
    """
    withoutdropping(ad::AdaptiveBase) = NoDropWrapper(ad)

    function _fit!(wrap::NoDropWrapper, value)
        ad = wrap.ad 
        fit!(ad.window[1], value)
        fit!(ad.stats, value)
        
        compress!(ad)   
        wrap
    end

    for fun in (:nobs, :value, :stats, :mean)
        @eval ($fun)(x::NoDropWrapper, args...) = ($fun(x.ad, args...))
    end

    struct MaxLength <: AdaptiveBase
        ad::AdaptiveBase
        maxlength::Int 
    end

    withmaxlength(ad::AdaptiveBase, maxlength) = MaxLength(ad, maxlength * M)

    function _fit!(wrap::MaxLength, value)
        _fit!(wrap.ad, value)
        if length(wrap.ad.window) > wrap.maxlength
            while length(wrap.ad.window) > wrap.maxlength
                pop!(wrap.ad.window)
            end
            wrap.ad.stats = Variance()
            for var in wrap.ad.window
                merge!(wrap.ad.stats, var)
            end
        end
    end

    for fun in (:nobs, :value, :stats, :mean)
        @eval ($fun)(x::MaxLength, args...) = ($fun(x.ad, args...))
    end


end # module

