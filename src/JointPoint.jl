# src/JointPoint.jl

module JointPoint

using Interpolations
using Optim

export findjoint

function findjoint(xc::AbstractVector{<:AbstractFloat}, 
				   yc::AbstractVector{<:AbstractFloat}, k::Integer; kwargs...)
	@assert issorted(xc, lt=<=)
	n = length(xc)
	@assert n == length(yc)
	xcl = xc[begin]
	xcr = xc[end]
	function p2xy(params)
		mat = [xcl; params[1:2*k+1]; xcr; params[2*k+2]]
		xj, yj = eachrow(sortslices(reshape(mat, 2, k+2), dims=2))
	end
	function f(params)
		xj, yj = p2xy(params)
		interp = linear_interpolation(xj, yj)
		return sqrt(sum((interp.(xc) .- yc) .^ 2))
	end
	ind = round.(Int, range(1, 10, length=k+2))
	params = permutedims(hcat(xc[ind], yc[ind]))[[2:2*k+2; 2*k+4]]
	niter = 3000
	while true
		result = optimize(f, params; iterations=niter, kwargs...)
		params = result.minimizer
		Optim.converged(result) && break
		niter *= 2
	end
	return collect.(p2xy(params))
end
findjoint(xc::AbstractVector{<:Real}, yc::AbstractVector{<:Real}, k::Integer) = 
	findjoint(float(xc), float(yc), k)

end # module JointPoint
