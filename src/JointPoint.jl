# src/JointPoint.jl

module JointPoint

using Interpolations
using LogExpFunctions
using Optim

export findjoint

function findjoint(xc::AbstractVector{<:AbstractFloat}, 
				   yc::AbstractVector{<:AbstractFloat}, k::Integer; kwargs...)
	@assert issorted(xc, lt=<=)
	n = length(xc)
	@assert n == length(yc)
	xcl = xc[begin]
	xcr = xc[end]
	u2x(u) = (xcr - xcl) * logistic(u) + xcl
	x2u(x) = logit((x - xcl) / (xcr - xcl))
	function p2xy(params)
		mat = [-Inf; params[1:2*k+1]; +Inf; params[2*k+2]]
		uj, yj = eachrow(sortslices(reshape(mat, 2, k+2), dims=2))
		xj = u2x.(uj)
		return xj, yj
	end
	function f(params)
		xj, yj = p2xy(params)
		interp = linear_interpolation(xj, yj)
		return sqrt(sum((interp.(xc) .- yc) .^ 2))
	end
	ind = round.(Int, range(1, 10, length=k+2))
	params = permutedims(hcat(x2u.(xc[ind]), yc[ind]))[[2:2*k+2; 2*k+4]]
	niter = 3000
	while true
		result = optimize(f, params; 
			x_tol=0.0, f_tol=0.0, g_tol=1e-12, iterations=niter, kwargs...)
		params = result.minimizer
		Optim.converged(result) && break
		niter *= 2
	end
	return collect.(p2xy(params))
end
findjoint(xc::AbstractVector{<:Real}, yc::AbstractVector{<:Real}, k::Integer) = 
	findjoint(float(xc), float(yc), k)

end # module JointPoint
