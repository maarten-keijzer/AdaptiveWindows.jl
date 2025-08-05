module AdaptiveWindows

export AdWin, AdWinGroup
export AdaptiveMean # backward compatibility
export fit!, value, mean, nobs, stats, withoutdropping, withmaxlength, update_no_check!

import StatsBase: nobs, fit!, merge!, var, mean
import OnlineStatsBase: value, OnlineStat, Variance, Mean, _fit!

#=
 Adaptive Windowing version 2 (AdaptiveMean2)

    Used to track the mean value of a stream of data with a possibly changing population.
    If the population is determined to be changed, older observations will be dropped

  Bifet and Gavalda. Learning from Time-Changing Data with Adaptive Windowing

  It is assumed that values are between 0 and 1 (TODO)
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
    mutable struct AdWin <: AdaptiveBase
        
        δ ::Float64
        window::Vector{Variance}
        stats ::Variance
        
        onshiftdetect

        AdWin(;δ = 0.001, onshiftdetected = noaction) = new(δ, [Variance() for _ in 1:M], Variance(), onshiftdetected)    
    end
    AdaptiveMean = AdWin

    noaction(ad, idx) = nothing

    function _fit!(ad::AdWin, value)
        fit!(ad.window[1], value)
        fit!(ad.stats, value)
    
        compress!(ad)
        idx = checkifdrifting(ad)
        if idx < typemax(Int)
            ad.onshiftdetect(ad, idx)
        end
        if idx < typemax(Int)
            drop!(ad, idx)
        end
        ad
    end

    function update_no_check!(ad::AdWin, value)
        fit!(ad.window[1], value)
        fit!(ad.stats, value)
        compress!(ad)
        ad
    end

    nobs(ad::AdWin) = ad.stats.n
    stats(ad::AdWin) = ad.stats

    function mean(ad::AdWin)
        ad.stats.μ
    end

    function value(ad::AdWin)
        mean(ad)
    end
                
    function compress!(ad::AdWin)
        makespace!(ad, 1, 1.0);
    end
    
    function print_trace(m)
        for v in m.window
            print(nobs(v), " ")
        end
        println("stats: ", nobs(m.stats))
    end    

    function makespace!(ad::AdWin, start::Int, max::Float64)
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

    function checkifdrifting(ad::AdWin)
        statsToRight = tomean(ad.stats)
        statsToLeft = Mean()
    
        deltaPrime = ad.δ / log(nobs(ad.stats))
        
        logDeltaPrime = log(2/deltaPrime)
        variance = var(ad.stats)
    

        for i in 1:length(ad.window)
            if nobs(ad.window[i]) == 0
                continue
            end
    
            _remove!(statsToRight, ad.window[i])
            _merge!(statsToLeft, ad.window[i])
            
            if statsToRight.n == 0
                break # only zeros from this point on
            end
    
            mInv = 1.0/ statsToRight.n + 1.0/ statsToLeft.n;
    
            # Less sensitive version:
            # epsCutOrg = sqrt(2 * mInv * log(4/deltaPrime))

            epsCut =sqrt(2 * mInv * variance * logDeltaPrime) + 2.0/3.0 * mInv * logDeltaPrime;
            
            if abs(statsToRight.μ - statsToLeft.μ) > epsCut
                # return the index of the first window to drop
                return i
            end
        end
    
        return typemax(Int)
    end
    
    function drop!(ad::AdWin, i::Int)
        # Drift detected, clear all from here on
        for j = i+1:length(ad.window)
            ad.window[j] = Variance()
        end
        ad.stats = Variance()
        for stat in 1:i
            merge!(ad.stats, ad.window[stat])
        end
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


    # struct AdaptiveMultinomial <: AdaptiveBase
    #     error_tracker::AdaptiveBase
    #     class_trackers::Vector{AdaptiveBase}
    # end

    # """
    #     Creates a set of AdWins that tracks the distribution of a class variable.
    #     The error tracker tracks the cross entropy of the class distribution.
    #     The class trackers track the mean of the class variable.
    #     The class trackers are updated with 1.0 if the class variable is the same as the class value,
    #     and 0.0 otherwise. The probabilities add up to 1.0.
    #     The error tracker is updated with the negative log of the mean of the class variable.

    #     When the distribution of the loss tracker changes, the class trackers are culled to the same length.
    # """
    # function AdaptiveMultinomial(nclasses::Int;δ = 0.001) 
    #     class_trackers = [withoutdropping(AdaptiveMean()) for _ in 1:nclasses]
    #     error_tracker = AdaptiveMean(δ = δ, onshiftdetected = sync_if_shifted(class_trackers))
    #     AdaptiveMultinomial(error_tracker, class_trackers)
    # end

    # function sync_if_shifted(class_trackers)
    #     (_, idx) -> begin 
    #         for class_tracker in class_trackers
    #             drop!(class_tracker.ad, idx)
    #         end
    #     end
    # end

    # function _fit!(mnAdWin::AdaptiveMultinomial, class_value::Int)
    #     @assert 0 < class_value <= length(mnAdWin.class_trackers)

    #     nats = 0.0;
    #     for i in eachindex(mnAdWin.class_trackers)
    #         if i == class_value
    #             mn = mean(mnAdWin.class_trackers[i])
    #             if 0.0 < mn < 1.0
    #                 nats -= log(mn)
    #             end
    #             fit!(mnAdWin.class_trackers[i], 1.0)
    #         else
    #             fit!(mnAdWin.class_trackers[i], 0.0)
    #         end
    #     end
    #     # track cross entropy of the class distribution
    #     fit!(mnAdWin.error_tracker, nats)
        
    # end

    # for fun in (:nobs, :value, :stats, :mean)
    #     @eval ($fun)(x::AdaptiveMultinomial, args...) = ($fun(x.error_tracker, args...))
    # end


    struct AdWinGroup <: AdaptiveWindows.AdaptiveBase 
        adwins::Vector{AdaptiveWindows.AdWin}
    end

    Base.getindex(syncedAdWins::AdWinGroup, i::Int) = syncedAdWins.adwins[i]

    function AdWinGroup(nadwins::Int;δ = 0.001, onshiftdetected = noaction) 
        adwins = [AdWin(δ = δ, onshiftdetected = onshiftdetected) for _ in 1:nadwins]
        AdWinGroup(adwins)
    end
    
    function fit!(syncedAdWins::AdWinGroup, values::Vector{T}; check_shift = true) where T <: Number
        @assert 0 < length(values) <= length(syncedAdWins.adwins)
    
        # fit each adwin, but don't drop just yet because windows might get out of sync
        for i in eachindex(syncedAdWins.adwins)
            adwin = syncedAdWins.adwins[i] 
            value = values[i]
            fit!(adwin.window[1], value)
            fit!(adwin.stats, value)
            
            compress!(adwin)   
        end
        
        if check_shift
        # drop to smallest window if any adwin detects a shift
            idx = typemax(Int)
            for ad in syncedAdWins.adwins
                idx = min(idx, checkifdrifting(ad))
            end

            if idx < typemax(Int)
                for ad in syncedAdWins.adwins
                    drop!(ad, idx)
                end
            end
        end
    end

    nobs(x::AdWinGroup) = nobs(x.adwins[1])
    for fun in (:value, :stats, :mean)
        @eval ($fun)(x::AdWinGroup) = ($fun.(x.adwins))
    end
    
end # module

