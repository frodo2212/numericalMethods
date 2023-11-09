include("../Matrixfunctions.jl")

#Die Matrix brauche ich eigentlich gar nicht in der Form - muss ja erst die ecken entfernt haben
matrix = zeros(Float64, 10, 10)
[matrix[i,i] = 2 for i in (1:10)]
[matrix[i,i-1] = 1 for i in (2:10)]
[matrix[i-1,i] = -1 for i in (2:10)]
matrix[10,1] = -1
matrix[1,10] = 1
b0 = ones(Float64, 10)
[b0[i] = (-1)^i*(i-1) for i in (1:10)]


function calculate_x(matrix::Matrix, vector::Vector)
    tridiagmatrix, u, v = matrix_tridag_zerlegen(matrix, 1)
    _,y = tridiag_x_calculation(tridiagmatrix, vector)
    _,z = tridiag_x_calculation(tridiagmatrix, u)
    x = y-(scalprod(v,y)/(1+scalprod(v,z)))*z
    if !isapprox(matrix*x, vector)
        error("solution not valid")
    end
    return x
end

function calculate_x_via_inversionmethod(matrix::Matrix, vector::Vector)
    tridiagmatrix, u, v = matrix_tridag_zerlegen(matrix, 1)
    inverse = tridiag_invert(tridiagmatrix)
    z = inverse*u
    w = inverse'*v
    lambda = scalprod(v,z)
    newinverse = inverse-(outprod(z,w))/(1+lambda)
    x = newinverse*vector
    if !isapprox(matrix*x, vector)
        error("solution not valid")
    end
    return x
end


x1 = calculate_x(matrix, b0)
matrix*x1
x2 = calculate_x_via_inversionmethod(matrix, b0)
matrix*x2


#Test über das allgemeine verfahren
inv = calculate_inverse(matrix)
