# test/runtests.jl

using Aqua
using JointPoint
using Test

Aqua.test_ambiguities(JointPoint)
Aqua.test_unbound_args(JointPoint)
Aqua.test_undefined_exports(JointPoint)
Aqua.test_piracies(JointPoint)
Aqua.test_project_extras(JointPoint)
Aqua.test_stale_deps(JointPoint)
Aqua.test_deps_compat(JointPoint)
Aqua.test_deps_compat(JointPoint)

@testset "findjoint" begin
	xc = 0:10
	yc = [0:2:10; 13:3:25]
	n = length(xc)
	a = (n * sum(xc .* yc) - sum(xc) * sum(yc)) / 
		(n * sum(xc .* xc) - sum(xc) ^ 2)
	b = n \ sum(yc) - a / n * sum(xc)
	xj0, yj0 = findjoint(xc, yc, 0)
	@test isapprox(xj0, [0, 10]; atol=1e-6)
	@test isapprox(yj0, a .* [0, 10] .+ b; atol=1e-3)
	xj1, yj1 = findjoint(xc, yc, 1)
	@test isa(xj1, Vector{Float64})
	@test isa(yj1, Vector{Float64})
	@test isapprox(xj1, [0, 5, 10]; atol=1e-6)
	@test isapprox(yj1, [0, 10, 25]; atol=1e-6)
end
