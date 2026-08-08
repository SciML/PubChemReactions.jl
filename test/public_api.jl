using Catalyst: Reaction, @species
using PubChemReactions
using Symbolics: @variables
using Test

compound_json = PubChemReactions.parse_json(
    """
    {
      "PC_Compounds": [{
        "id": {"id": {"cid": 962}},
        "atoms": {"aid": [1, 2, 3], "element": [8, 1, 1]},
        "bonds": {"aid1": [1, 1], "aid2": [2, 3]},
        "charge": 0
      }]
    }
    """
)
view_json = PubChemReactions.parse_json(
    """
    {
      "Record": {
        "RecordTitle": "Water",
        "RecordNumber": 962,
        "Section": []
      }
    }
    """
)

@variables t
water = only(@species water(t))
water = PubChemReactions.tospecies(water; jsons = (compound_json, view_json))

@test isspecies(water)
@test get_cid(water) == 962
@test get_name(water) == "Water"
@test get_charge(water) == 0
@test atom_counts(water) == Dict(8 => 1, 1 => 2)
@test sum(values(element_counts(water))) == 3
@test get_graph(water) !== nothing
@test occursin("cid = 962", make_at_species(water))

reaction = Reaction(1, [water], [water])
@test isbalanced(reaction)
@test length(balance_eqs(reaction)) == 3
@test size(atom_matrix(reaction), 2) == 2
@test eqs_to_mathematica(balance_eqs(reaction)) isa String

mktempdir() do path
    saved = save_species(water; path)
    @test saved == joinpath(path, "962")
    @test isfile(joinpath(saved, "pug.json"))
    @test isfile(joinpath(saved, "pug_view.json"))
end
