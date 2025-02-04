a = [1 2 3; 4 5 6; 7 8 9; 10 11 15]

for idx in CartesianIndices(a)
    println(idx)
    i, j = Tuple(idx)  # Extract row and column indices
    println("a[$i, $j] = $(a[i, j])")
end



model = ExaCore(Float64; backend = nothing)
x = variable(model, 4, 3)
obj = objective(model, x[i, j] for i in 1:4, j in 1:3)
c1 = constraint(model, x[i, j] for i in 1:4, j in 1:3; lcon = a, ucon = fill(Inf, 3, 3))

result = ipopt(ExaModel(model))
println(solution(result, x))