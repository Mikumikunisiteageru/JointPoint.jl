# src/ZigZagFit.jl

module ZigZagFit

using Interpolations
using Optim

export zigzagfit

function zigzagfit(xx0::AbstractVector{<:AbstractFloat}, 
				   yy0::AbstractVector{<:AbstractFloat}, k::Integer; kwargs...)
	@assert issorted(xx0, lt=<=)
	n = length(xx0)
	@assert n == length(yy0)
	xl = xx0[begin]
	xr = xx0[end]
	function f(xyparams)
		xyparammat = [xl; xyparams[1:2*k+1]; xr; xyparams[2*k+2]]
		xxp, yyp = eachrow(sortslices(reshape(xyparammat, 2, k+2), dims=2))
		interp = linear_interpolation(xxp, yyp)
		return sqrt(sum((interp.(xx0) .- yy0) .^ 2))
	end
	ii = round.(Int, range(1, 10, length=k+2))
	xxt = xx0[ii]
	yyt = yy0[ii]
	xyparamst = permutedims([xxt yyt])[[2:2*k+2; 2*k+4]] 
	result = optimize(f, xyparamst; iterations=30000, kwargs...)
	@assert Optim.converged(result)
	xyparamsm = result.minimizer
	xyparammatm = [xl; xyparamsm[1:2*k+1]; xr; xyparamsm[2*k+2]]
	xxm, yym = collect.(eachrow(sortslices(reshape(xyparammatm, 2, k+2), dims=2)))
	return xxm, yym
end
zigzagfit(xx0::AbstractVector{<:Real}, yy0::AbstractVector{<:Real}, k::Integer) = 
	zigzagfit(float(xx0), float(yy0), k)

end # module ZigZagFit
