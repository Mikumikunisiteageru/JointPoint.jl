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

struct IntervalSum{T<:AbstractFloat}
	c::Vector{T}
	IntervalSum(v::AbstractVector{T}) where {T <: AbstractFloat} = new{T}(cumsum(v))
end

Base.getindex(f::IntervalSum, a::Integer, b::Integer) = a > 0 ? f.c[b] - f.c[a] : f.c[b]

function findjoint(xc::AbstractVector{T}, 
				   yc::AbstractVector{T}, k::Integer; kwargs...) where {T<:AbstractFloat}
	issorted(xc, lt=<=) || throw(ArgumentError("`xc` must be strictly increasing"))
	n = length(xc)
	n == length(yc) || throw(ArgumentError("`xc` and `yc` must have the same length"))
	cxi2  = IntervalSum(xc .^ 2)
	cxi   = IntervalSum(xc)
	cxiyi = IntervalSum(xc .* yc)
	cyi   = IntervalSum(yc)
	cone  = IntervalSum(ones(T, n))
	m = k + 2
	function xjcore2xjyj(xjcore)
		xj = [xc[begin]; sort(xjcore); xc[end]]
		dl = zeros(T, m-1)
		d  = zeros(T, m)
		du = zeros(T, m-1)
		right = zeros(T, m)
		# creating new vectors are faster than zeroing old ones since they are small
		id = [0; searchsortedlast.([xc], xj[2:m])]
		for j = 1:m-1
			sm1, t = id[j], id[j+1]
			sxi2  = cxi2[sm1, t]
			sxi   = cxi[sm1, t]
			sxiyi = cxiyi[sm1, t]
			syi   = cyi[sm1, t]
			sone  = cone[sm1, t]
			a, b = xj[j:j+1]
			bma = b - a
			d[j:j+1]  .+= (a+b) * sxi - sxi2 - a*b * sone
			du[j]       = sxi2 - 2*a * sxi + a^2 * sone
			dl[j]       = sxi2 - 2*b * sxi + b^2 * sone
			right[j]   += bma * (sxiyi - a * syi)
			right[j+1] += bma * (b * syi - sxiyi)
		end
		left = Tridiagonal(dl, d, du)
		yj = try
			left \ right
		catch
			pinv(left) * right
		end
		return xj, yj
	end
	function f(xjcore)
		xj, yj = xjcore2xjyj(xjcore)
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
	return xjcore2xjyj(xjcore)
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
