# matlab
u = symmatrix('u', [3 1])
x = symmatrix2sym(u)
eqs = [
    x(1) == x(3),
    x(1) + x(2) == 2 * x(3),
    -x(1) + x(2) == 0,
    x(1) > 0
]
solution = solve(eqs,[sym('u1'),sym('u2'),sym('u3')],...
    'PrincipalValue',true,...
    'ReturnConditions',true);

solution.u1, solution.u2, solution.u3 % (1,1,1)

# matlab 2
```Matlab
syms x [3 1] positive integer
eqs = [
    x(1) == x(3)
    x(1) + x(2) == 2 * x(3)
    -x(1) + x(2) == 0
]

solution = solve(eqs,[sym('x1'),sym('x2'),sym('x3')])
[solution.x1, solution.x2, solution.x3] % [1,1,1]
```

# mathematica
eqs = {
    u[1] == u[3],
    u[1] + u[2] == 2*u[3],
    -u[1] + u[2] == 0
}
possible = FindInstance[eqs, {u[1], u[2], u[3]}, PositiveIntegers]
xs = Map[Last, First[possible]]
sol = xs / Min[xs] # {1,1,1}

# python 
eqs = [
    Eq(x, z),
    Eq(x + y, 2 * z),
    Eq(-x + y, 0)
]
solve(eqs, (x, y, z), domain=S.Naturals)



# julia 
using Symbolics
@variables u[1:3]
eqs = [
    u[1] ~ u[3],
    u[1] + u[2] ~ 2 * u[3],
    -u[1] + u[2] ~ 0
]



ts = PubChemReactions.eq_to_term.(eqs)
a, b, isl = linear_expansion(ts, collect(x[1:4]))
a
h, w = size(a)
A = zeros(h, w + 1)
A[1:h, 1:w] .= value.(a)
@test isbalanced(rxn)
@test_throws balance(rxn)
@variables u[1:3]

eqs3 = [
    u[1] ~ u[3],
    u[1] + u[2] ~ 2 * u[3],
    -u[1] + u[2] ~ 0
]

seqs = substitute(eqs3, Dict([u[3] => 1]))
Symbolics.solve_for(seqs, u; simplify=false)


A = [
    1 0 -1
    1 1 -2
    -1 1 0
]
collect(A * u)


A = [
    1 0 -1 0
    1 1 -2 0
    -1 1 0 0
    0 0 0 1
]
Symbolics.solve_for(eqs3, u; simplify=false)
eqs4 = [
    u[1] + u[2] ~ 2
    u[2] - u[1] ~ 0
]
Symbolics.solve_for(eqs4, u[1:2]; simplify=false)


eqs = [u[1] ~ u[2], u[2] ~ u[1]]
solve_for(eqs, var)

s = PubChemReactions.parse_rhea_equation("")

eqs = balance_eqs(rxn)
ts = eq_to_term.(eqs)
vars = unique(reduce(vcat, Symbolics.get_variables.(ts))) # get_variables permutes, since the order they show in eqs
# vars = unique(reduce(vcat, x[1:4]))
vars = collect(x[1:4])
a, b, islinear = linear_expansion(eq_to_term.(balance_eqs(rxn)), vars)

function eq_str_to_wl(str)
    str = replace(str, "~" => "==")
    # "{$str}"
end

function eqs_to_mathematica(eqs)
    es = string.(eqs)
    eq_strs = map(eq_str_to_wl, es)
    join(["{", join(eq_strs, ", "), "}"])
end
    