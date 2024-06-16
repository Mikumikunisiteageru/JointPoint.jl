# src/JointPoint.jl

module JointPoint

using BlackBoxOptim
using Interpolations

export findjoint, aicc

function rss(xc::AbstractVector{<:AbstractFloat}, yc::AbstractVector{<:AbstractFloat}, 
			 xj::AbstractVector{<:AbstractFloat}, yj::AbstractVector{<:AbstractFloat})
	interp = linear_interpolation(xj, yj)
	return sum((interp.(xc) .- yc) .^ 2)
end

function findjoint(xc::AbstractVector{<:AbstractFloat}, 
				   yc::AbstractVector{<:AbstractFloat}, k::Integer; p=1, kwargs...)
	@assert issorted(xc, lt=<=)
	n = length(xc)
	@assert n == length(yc)
	xcl = xc[begin]
	xcr = xc[end]
	ycd, ycu = extrema(yc)
	function p2xy(params)
		mat = [xcl; params[1:2*k+1]; xcr; params[2*k+2]]
		return eachrow(sortslices(reshape(mat, 2, k+2), dims=2))
	end
	function f(params)
		xj, yj = p2xy(params)
		return rss(xc, yc, xj, yj)
	end
	ind = round.(Int, range(1, 10, length=k+2))
	params = permutedims(hcat(xc[ind], yc[ind]))[[2:2*k+2; 2*k+4]]
	xrange = (xcl, xcr)
	yrange = ((p+1) * ycd - p * ycu, (p+1) * ycu - p * ycd)
	searchrange = vcat(fill(xrange, 1, k+2), fill(yrange, 1, k+2))[[2:2*k+2; 2*k+4]]
	res = bboptimize(f; SearchRange=searchrange, MaxStep=300000, MaxTime=10, 
		MinDeltaFitnessTolerance=1e-20, kwargs...)
	params = best_candidate(res)
	return collect.(p2xy(params))
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
