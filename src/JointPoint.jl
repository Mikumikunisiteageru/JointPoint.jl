# src/JointPoint.jl

module JointPoint

using BlackBoxOptim
using Interpolations

export findjoint

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
		interp = linear_interpolation(xj, yj)
		return sqrt(sum((interp.(xc) .- yc) .^ 2))
	end
	ind = round.(Int, range(1, 10, length=k+2))
	params = permutedims(hcat(xc[ind], yc[ind]))[[2:2*k+2; 2*k+4]]
	xrange = (xcl, xcr)
	yrange = ((p+1) * ycd - p * ycu, (p+1) * ycu - p * ycd)
	searchrange = vcat(fill(xrange, 1, k+2), fill(yrange, 1, k+2))[[2:2*k+2; 2*k+4]]
	res = bboptimize(f; SearchRange=searchrange, TraceMode=:silent, kwargs...)
	params = best_candidate(res)
	return collect.(p2xy(params))
end
findjoint(xc::AbstractVector{<:Real}, yc::AbstractVector{<:Real}, k::Integer) = 
	findjoint(float(xc), float(yc), k)

end # module JointPoint
