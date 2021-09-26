C6H12O6 = PubChemReactions.gen_sym("glucose")
H2O = PubChemReactions.gen_sym("water")

res = PubChemReactions.species_info(C6H12O6)
res[:InformationList][:Information]

df = PubChemReactions.get_biochem_rxns(C6H12O6, H2O)
eqs = df[!,:Equation]
eq = first(eqs) # I happen to know this one has no prefixed stoich data
reactants, products, rstoich, pstoich = PubChemReactions.parse_rhea_equation(eq)

# arbitrarily assigning constant rate of 1. 
# this is where we'd like to do some lookup on the reactants and products to make an educated guess about the rate law
rxn = Reaction(1, reactants, products, rstoich, pstoich; only_use_rate=true) 

@test rxn isa Catalyst.Reaction

eq = eqs[19] # includes stoich 
rxn = Reaction(1, PubChemReactions.parse_rhea_equation(eq)...; only_use_rate=true) 
@test rxn isa Catalyst.Reaction

eqs = eqs[1:20]
rxns = Reaction[]
@sync @async for (i, eq) in enumerate(eqs)
    @info i, eq
    try
        res = PubChemReactions.parse_rhea_equation(eq)
        res === nothing && continue
        rxn = Reaction(1, res...; only_use_rate=true) 
        push!(rxns, rxn)
        @info rxn
    catch e 
        @info e
    end
end

N = length(rxns)
@info "$N / $(length(eqs))"
# @test N == 57

h2o = rxns[1].substrates[1]
h2o_ = rxns[2].substrates[2]
@test isequal(h2o, h2o_)


sys_states = rxns_to_states(rxns)
@named mysys = ReactionSystem(rxns, Catalyst.DEFAULT_IV, sys_states, []; defaults=Dict(sys_states .=> 1))
sys = convert(ODESystem, mysys)

# equations(sys)
# prob = ODEProblem(sys, [], (0, 100.))
# sol = solve(prob, Tsit5())
