# src/JointPoint.jl

module JointPoint

using BlackBoxOptim
using Interpolations
using LinearAlgebra

export findjoint, aicc

function rss(xc, yc, xj, yj)
	interp = linear_interpolation(xj, yj)
	return sum((interp.(xc) .- yc) .^ 2)
end

function xjcore2xjyj(xc, yc, xjcore)
	T = float(promote_type(eltype(xc), eltype(yc)))
	xj = [xc[begin]; xjcore; xc[end]]
	m = length(xj)
	dl = zeros(T, m-1)
	d = zeros(T, m)
	du = zeros(T, m-1)
	right = zeros(T, m)
	id = [0; searchsortedlast.([xc], xj[2:m])]
	for j = 1:m-1
		i = id[j]+1:id[j+1]
		a, b = xj[j:j+1]
		bma = b - a
		bmx = b .- xc[i]
		xma = xc[i] .- a
		d[j:j+1] .+= sum(bmx .* xma)
		du[j] = sum(xma .^ 2)
		dl[j] = sum(bmx .^ 2)
		right[j] += bma * sum(xma .* yc[i])
		right[j+1] += bma * sum(bmx .* yc[i])
	end
	left = Tridiagonal(dl, d, du)
	yj = try
		left \ right
	catch
		pinv(left) * right
	end
	return xj, yj
end

function findjoint(xc::AbstractVector{T}, 
				   yc::AbstractVector{T}, k::Integer; kwargs...) where {T<:AbstractFloat}
	issorted(xc, lt=<=) || throw(ArgumentError("`xc` must be strictly increasing"))
	n = length(xc)
	n == length(yc) || throw(ArgumentError("`xc` and `yc` must have the same length"))
	function f(xjcore)
		xj, yj = xjcore2xjyj(xc, yc, sort(xjcore))
		return rss(xc, yc, xj, yj)
	end
	xjcore = if k > 0
		searchrange = fill((xc[begin], xc[end]), k)
		res = bboptimize(f; SearchRange=searchrange, TraceMode=:silent, MaxTime=5*k, kwargs...)
		startswith(BlackBoxOptim.stop_reason(res), "Max time") && 
			@warn("optimization process may not have converged")
		best_candidate(res)
	else
		T[]
	end
	return xjcore2xjyj(xc, yc, sort(xjcore))
end
findjoint(xc::AbstractVector{<:AbstractFloat}, 
		  yc::AbstractVector{<:AbstractFloat}, k::Integer; kwargs...) = 
	findjoint(promote(xc, yc)..., k; kwargs...)
findjoint(xc::AbstractVector{<:Real}, yc::AbstractVector{<:Real}, k::Integer; kwargs...) = 
	findjoint(float(xc), float(yc), k; kwargs...)

function aicc(xc::AbstractVector{T}, yc::AbstractVector{T}, 
			  xj::AbstractVector{T}, yj::AbstractVector{T}) where {T<:AbstractFloat}
	n = length(xc)
	karg = length(xj) * 2 - 2
	return 2 * karg * n / (n - karg - 1) + n * log(rss(xc, yc, xj, yj) / n)
end
aicc(xc::AbstractVector{<:AbstractFloat}, yc::AbstractVector{<:AbstractFloat}, 
	 xj::AbstractVector{<:AbstractFloat}, yj::AbstractVector{<:AbstractFloat}) = 
	 	aicc(promote(xc, yc, xj, yj)...)
aicc(xc::AbstractVector{<:Real}, yc::AbstractVector{<:Real}, 
	 xj::AbstractVector{<:Real}, yj::AbstractVector{<:Real}) = 
		aicc(float(xc), float(yc), float(xj), float(yj))

end # module JointPoint
