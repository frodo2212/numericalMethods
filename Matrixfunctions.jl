function scalprod(a::Vector, b::Vector)
    if length(a)!=length(b)
        error("dimensions dont match")
    end
    out =0
    [out+=(a[i]*b[i]) for i in (1:length(a))]
    return out
end

function outprod(a::Vector, b::Vector)
    len = length(a)
    if len!=length(b)
        error("dimensions dont match")
    end
    out = zeros(Float64, len, len)
    [out[i,j] = a[i]*b[j] for i in (1:len), j in (1: len)]
    return out
end

function backward_substitution(input::Matrix, vector::Vector) #passt so
    upper_matrix  = deepcopy(input)
    vector2 = deepcopy(vector)
    #ab hier die substitution
    rows, collumns = size(upper_matrix)
    out = Vector{Float64}(undef, rows)
    for x in round.(Int64,collect(LinRange(rows,1, rows)))
        subtraction = 0
        for i in (x+1:rows)
            subtraction += upper_matrix[x,i]*out[i] 
        end
        out[x] = (vector2[x]-subtraction)/upper_matrix[x,x]
    end
    return out
end

function forward_substitution(input::Matrix, vector::Vector)
    upper_matrix  = deepcopy(input)
    vector2 = deepcopy(vector)
    rows, collumns = size(upper_matrix)
    out = Vector{Float64}(undef, rows)
    for x in (1:rows)
        subtraction = 0
        for i in (1:x-1)
            subtraction += upper_matrix[x,i]*out[i] 
        end
        out[x] = vector2[x]-subtraction
    end
    return out
end

function GaußJordan_wo_pivot(input::Matrix{Float64}, vec_input::Vector{Float64})
    vector = deepcopy(vec_input)
    upper_matrix = deepcopy(input)
    rows, collumns = size(inv_matrix)
    if rows != collumns
        error("Matrix is not squared")
    end
    for a in (1:rows-1)
        factor = upper_matrix[a,a] 
        for i in (1:rows)
            upper_matrix[a,i] = upper_matrix[a,i] / factor
        end
        vector[a] = vector[a] / factor
        for c in (a+1:rows)            
            subtraction_factor = upper_matrix[c,a]
            for i in (1:rows)
                upper_matrix[c,i] = upper_matrix[c,i]-subtraction_factor*upper_matrix[a,i] 
            end
            vector[c] = vector[c]-subtraction_factor*vector[a]
        end
    end
    output = backward_substitution(upper_matrix, vector)
    if !isapprox(input*output, vec_input)
        error("calculation without pivoting didn't work")
    end
    return output
end

function GaußJordan_partial_pivot(input::Matrix{Float64}, vec_input::Vector{Float64})
    vector = deepcopy(vec_input)
    upper_matrix = deepcopy(input)
    rows, collumns = size(upper_matrix)
    if rows != collumns
        error("Matrix is not squared")
    end
    for a in (1:rows-1)
        max = argmax(abs.(upper_matrix[a:rows, a]))+a-1
        tmpvec = upper_matrix[a,:]
        upper_matrix[a,:] = upper_matrix[max,:] 
        upper_matrix[max,:] = tmpvec
        tmpval = vector[a]
        vector[a] = vector[max]
        vector[max] = tmpval
        factor = upper_matrix[a,a] 
        for i in (1:rows)
            upper_matrix[a,i] = upper_matrix[a,i] / factor
        end
        vector[a] = vector[a] / factor
        for c in (a+1:rows)            
            subtraction_factor = upper_matrix[c,a]
            for i in (1:rows)
                upper_matrix[c,i] = upper_matrix[c,i]-subtraction_factor*upper_matrix[a,i] 
            end
            vector[c] = vector[c]-subtraction_factor*vector[a]
        end
    end
    output = backward_substitution(upper_matrix, vector)
    if !isapprox(input*output, vec_input)
        error("calculation with pivoting didn't work")
    end
    return output
end

function calculate_inverse(A::Matrix) #TODO: das geht effizienter, wenn ich das beim Gauß direkt mit mache
    rows, collumns = size(A)
    if rows != collumns
        error("Matix not squared")
    end
    A_inverse = Matrix{Float64}(undef, rows, collumns)
    for r in (1:rows)
        b = Float64[i==r for i in (1:rows)]
        A_inverse[:,r] = GaußJordan_partial_pivot(A, b)
    end
    unitmatrix = zeros(Float64, rows, collumns)
    [unitmatrix[i,i] = 1 for i in (1:rows)]
    if !isapprox(A_inverse*A,unitmatrix)
        error("The solution is not a inverse")
    end
    return A_inverse
end

function crout(input::Matrix) 
    rows, collumns = size(input)
    out=Matrix{Float64}(undef, rows, collumns)
    for j in (1:rows)
        for i in (1:j)
            sum = 0
            for k in (1:i-1)
                sum += out[i,k]*out[k,j] 
            end
            out[i,j] = input[i,j] - sum
        end
        for i in (j+1:rows)
            sum = 0
            for k in (1:j)
                sum += out[i,k]*out[k,j] 
            end
            out[i,j] = (1/out[j,j])*(input[i,j]-sum)
        end
    end
    return out
end


#Für tridiag matrizen mit cyclic Randbedingungen
function matrix_tridag_zerlegen(matrix::Matrix, gamma = 0)
    if gamma == 0
        gamma = -matrix[1,1]
    end
    rows, collumns = size(matrix)
    u = zeros(Float64, rows)
    u[1] = gamma
    u[rows] = matrix[10,1]
    v = zeros(Float64, rows)
    v[1] = 1
    v[rows] = matrix[1,10]/gamma
    tridiagMatrix = deepcopy(matrix)
    tridiagMatrix[1,1] = matrix[1,1]-gamma
    tridiagMatrix[10,10] = matrix[10,10]-(matrix[10,1]*matrix[1,10])/gamma
    tridiagMatrix[1,10] = 0
    tridiagMatrix[10,1] = 0
    return tridiagMatrix, u, v
end

#TODO: besseren Namen überlegen
function tridiag_x_calculation(input::Matrix, vector::Vector)
    rows, collumns = size(input)
    y = Vector{Float64}(undef, rows)
    x = Vector{Float64}(undef, rows)
    out = Matrix{Float64}(undef, rows, collumns)
    y[1] = vector[1]
    out[1,1] = input[1,1]
    for i in (2:rows)
        out[i-1,i]=input[i-1,i]
        out[i,i-1] = input[i,i-1]/out[i-1,i-1]
        out[i,i] = input[i,i]-out[i,i-1]*out[i-1,i]
        y[i] = vector[i]-out[i,i-1]*y[i-1]
    end
    x[rows] = y[rows]/out[rows,rows]
    for i in round.(Int64,collect(LinRange(rows-1,1, rows-1)))
        x[i] = (y[i]-(out[i,i+1]*x[i+1]))/out[i,i]
    end
    return out, x
end

function tridiag_invert(input::Matrix)#über GaußJordan
    matrix = deepcopy(input)
    rows, collumns = size(input)
    out = zeros(Float64, rows, collumns)
    [out[i,i] = 1 for i in (1:rows)]
    for a in (1:rows)
        factor = matrix[a,a] 
        for i in (1:rows)
            matrix[a,i] = matrix[a,i] / factor
            out[a,i] = out[a,i] / factor
        end     
        if a < rows      
            subtraction_factor = matrix[a+1,a]
            for i in (1:rows)
                matrix[a+1,i] = matrix[a+1,i]-subtraction_factor*matrix[a,i] 
                out[a+1,i] = out[a+1,i]-subtraction_factor*out[a,i] 
            end
        end
    end
    for a in round.(Int64,collect(LinRange(rows,2, rows-1)))   
        subtraction_factor = matrix[a-1,a]
        for i in (1:rows)
            matrix[a-1,i] = matrix[a-1,i]-subtraction_factor*matrix[a,i] 
            out[a-1,i] = out[a-1,i]-subtraction_factor*out[a,i] 
        end
    end
    unit = zeros(Float64, rows, collumns)
    [unit[i,i] = 1 for i in (1:rows)]
    if !isapprox(matrix, unit)
        error("something went wrong")
    end
    return out
end