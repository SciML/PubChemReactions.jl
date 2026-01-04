using PubChemReactions, Catalyst, OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks

# str = "H2O2 = H2O + O2"
# s, p = PubChemReactions.parse_pathway_reaction(str)
@parameters t T
@variables A(t) B(t) C(t) D(t)
@parameters k1 = 10 k2 = 1
rxn = Reaction(k1, [A, B], [C, D])
rxn2 = Reaction(k2, [C, D], [A, B])
@named rs = ReactionSystem([rxn, rxn2], t)
sys = convert(ODESystem, rs)
u0 = [1, 1, 1, 1]
prob = ODEProblem(sys, u0, (0, 1.0e5))
cb = TerminateSteadyState(1.0e-4, 1.0e-3)
sol = OrdinaryDiffEq.solve(prob, Tsit5(); callback = cb)
t1 = sol.t[end]

prob2 = ODEProblem(sys, u0, (0, 1.0e5), (100, 10))
sol2 = OrdinaryDiffEq.solve(prob2, Tsit5(); callback = cb)
t2 = sol2.t[end]

# the equilibrium is the same
@test isapprox(sol[end], sol2[end]; atol = 1.0e-4)
@test t1 > 5 * t2
