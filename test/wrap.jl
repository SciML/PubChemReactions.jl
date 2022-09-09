

# struct Species end

# Symbolics.option_to_metadata_type(::Val{:species}) = Species
# source = getmetadata(Water, Symbolics.VariableSource)[1]
# T = Symbolics.option_to_metadata_type(Val(source))

# # Num{Species}(x::Num) = 
# abstract type Num2 <: Real end
# struct Num2{T} <: Num
#     val
# end
# ex = :(
#     struct Num2{T} <: Num
#         val
#     end
# )

# ex2 = :(struct Foo end)

# Species2 = Tuple{Num}