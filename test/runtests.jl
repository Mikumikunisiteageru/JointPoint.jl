# test/runtests.jl

using Aqua
using ZigZagFit
using Test

Aqua.test_ambiguities(ZigZagFit)
Aqua.test_unbound_args(ZigZagFit)
Aqua.test_undefined_exports(ZigZagFit)
Aqua.test_piracies(ZigZagFit)
Aqua.test_project_extras(ZigZagFit)
Aqua.test_stale_deps(ZigZagFit)
Aqua.test_deps_compat(ZigZagFit)
Aqua.test_deps_compat(ZigZagFit)

@testset "zigzagfit" begin
	xx0 = 0:10
	yy0 = [0:2:10; 13:3:25]
	n = length(xx0)
	a = (n * sum(xx0 .* yy0) - sum(xx0) * sum(yy0)) / 
		(n * sum(xx0 .* xx0) - sum(xx0) ^ 2)
	b = n \ sum(yy0) - a / n * sum(xx0)
	xxm0_ = [0, 10]
	yym0_ = a .* xxm0_ .+ b
	xxm0, yym0 = zigzagfit(xx0, yy0, 0)
	@test isapprox(xxm0, xxm0_; atol=1e-6)
	@test isapprox(yym0, yym0_; atol=1e-3)
	xxm1, yym1 = zigzagfit(xx0, yy0, 1)
	@test isa(xxm1, Vector{Float64})
	@test isa(yym1, Vector{Float64})
	@test isapprox(xxm1, [0, 5, 10]; atol=1e-6)
	@test isapprox(yym1, [0, 10, 25]; atol=1e-6)
end
