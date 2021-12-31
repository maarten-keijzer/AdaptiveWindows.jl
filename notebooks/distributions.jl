### A Pluto.jl notebook ###
# v0.17.4

using Markdown
using InteractiveUtils

# ╔═╡ 2a4e15a3-9b6c-42ce-9555-6a2e51d421e5
import Pkg; Pkg.add(url="https://github.com/maarten-keijzer/AdaptiveWindow.jl")

# ╔═╡ d22707fc-6a27-11ec-3e85-772cedeb1c06
begin
	using OnlineStats
	using AdaptiveWindow
	using Plots
end

# ╔═╡ 47278b8c-c601-4fdf-a687-06349db8f9aa
begin
	m1 = AdaptiveMean(δ = 0.001)
	m2 = Mean()
	m3 = Mean(weight=ExponentialWeight(.01))
	means1 = []
	means2 = []
	means3 = []
	for i = 1:1000
		v = i < 500 ? randn() : randn() + 1
		fit!(m1, v)
		fit!(m2, v)
		fit!(m3, v)
		push!(means1, value(m1))
		push!(means2, value(m2))
		push!(means3, value(m3))
	end
end

# ╔═╡ 4b10f740-dfb6-4e61-915d-e2d3226fd394
plot([means1, means2, means3], 
	label=["Adaptive Mean" "Regular Mean" "Exponential Moving Average"], 		 
    title="Distribution shifts at t=500", legend=:topleft)

# ╔═╡ Cell order:
# ╠═2a4e15a3-9b6c-42ce-9555-6a2e51d421e5
# ╠═d22707fc-6a27-11ec-3e85-772cedeb1c06
# ╠═47278b8c-c601-4fdf-a687-06349db8f9aa
# ╠═4b10f740-dfb6-4e61-915d-e2d3226fd394
