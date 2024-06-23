# src/JointPoint.jl

module JointPoint

using BlackBoxOptim
using Interpolations
using LinearAlgebra

export findjoint, aicc

function rss(xc::AbstractVector{<:AbstractFloat}, yc::AbstractVector{<:AbstractFloat}, 
			 xj::AbstractVector{<:AbstractFloat}, yj::AbstractVector{<:AbstractFloat})
	interp = linear_interpolation(xj, yj)
	return sum((interp.(xc) .- yc) .^ 2)
end

function xjinner2xjyj(xc, yc, xjinner)
	T = float(promote_type(eltype(xc), eltype(yc)))
	xj = [xc[begin]; xjinner; xc[end]]
	m = length(xj)
	dl = zeros(T, m-1)
	d = zeros(T, m)
	du = zeros(T, m-1)
	right = zeros(T, m)
	id = [0; searchsortedlast.([xc], xj[2:m])]
	for j = 1:m-1
		i = id[j]+1:id[j+1]
		a, b = xj[j:j+1]
		d[j:j+1] .+= sum(@.((b-xc[i]) * (xc[i]-a)))
		du[j] = sum(@.((xc[i]-a) ^ 2))
		dl[j] = sum(@.((b-xc[i]) ^ 2))
		right[j] += sum(@.((xc[i]-a) * (b-a) * yc[i]))
		right[j+1] += sum(@.((b-xc[i]) * (b-a) * yc[i]))
	end
	left = Tridiagonal(dl, d, du)
	yj = pinv(left) * right
	return xj, yj
end

function findjoint(xc::AbstractVector{<:AbstractFloat}, 
				   yc::AbstractVector{<:AbstractFloat}, k::Integer; kwargs...)
	@assert issorted(xc, lt=<=)
	n = length(xc)
	@assert n == length(yc)
	function f(xjinner)
		xj, yj = xjinner2xjyj(xc, yc, sort(xjinner))
		return rss(xc, yc, xj, yj)
	end
	xjinner = if k > 0
		searchrange = fill((xc[begin], xc[end]), k)
		res = bboptimize(f; SearchRange=searchrange, MaxStep=300000, MaxTime=10, 
			MinDeltaFitnessTolerance=1e-20, kwargs...)
		best_candidate(res)
	else
		eltype(xc)[]
	end
	return xjinner2xjyj(xc, yc, sort(xjinner))
end
findjoint(xc::AbstractVector{<:Real}, yc::AbstractVector{<:Real}, k::Integer) = 
	findjoint(float(xc), float(yc), k)

function aicc(xc::AbstractVector{<:AbstractFloat}, yc::AbstractVector{<:AbstractFloat}, 
			  xj::AbstractVector{<:AbstractFloat}, yj::AbstractVector{<:AbstractFloat})
	n = length(xc)
	karg = length(xj) * 2 - 2
	return 2 * karg * n / (n - karg - 1) + n * log(rss(xc, yc, xj, yj) / n)
end
aicc(xc::AbstractVector{<:Real}, yc::AbstractVector{<:Real}, 
	 xj::AbstractVector{<:Real}, yj::AbstractVector{<:Real}) = 
		aicc(float(xc), float(yc), float(xj), float(yj))

end # module JointPoint
