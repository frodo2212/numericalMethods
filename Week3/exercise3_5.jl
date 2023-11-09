matrix = Float64[4 1 5 3 7 3 5 1;12 31 1 5 7 2 17 3;1 2 3 4 5 6 7 4;6 3 7 11 10 9 3 7;3 9 8 3 2 1 4 1;4 15 12 3 17 1 1 4;
                5 2 12 14 12 12 13 1;9 9 7 5 3 1 1 2]
b1 = [1,4,5,8,9,3,1,3]
b2 = [5,0,1,3,7,0,0,7]
b3 = [1,0,1,0,1,0,6,5]

include("../Matrixfunctions.jl")

function do_all_the_work(matrix::Matrix, vectors::Vector{Vector{T}}) where T <: Real
    test = crout(matrix)
    output = Vector{Float64}[]
    for vec in vectors
        x = backward_substitution(test, forward_substitution(test, vec))
        push!(output, x)
        if !isapprox(matrix*x,vec)
            println("something went wrong")
        end
    end
    return output
end

y1, y2, y3 = do_all_the_work(matrix, [b1,b2,b3])
y3


matrix*y1
test = crout(matrix)
