function GaußJordan_wo_pivot(input::Matrix{Float64}, vec_input::Vector{Float64}) #funktioniert auch, wenn kein a,a wert null ist
    vector = deepcopy(vec_input)
    inv_matrix = deepcopy(input)
    rows, collumns = size(inv_matrix)
    if rows != collumns
        error("Matrix is not squared")
    end
    for a in (1:rows-1)
        factor = inv_matrix[a,a] 
        for i in (1:rows)
            inv_matrix[a,i] = inv_matrix[a,i] / factor
        end
        vector[a] = vector[a] / factor
        for c in (a+1:rows)            
            subtraction_factor = inv_matrix[c,a]
            for i in (1:rows)
                inv_matrix[c,i] = inv_matrix[c,i]-subtraction_factor*inv_matrix[a,i] 
            end
            vector[c] = vector[c]-subtraction_factor*vector[a]
        end
    end
    return inv_matrix, vector
end

function GaußJordan_w_pivot(input::Matrix{Float64}, vec_input::Vector{Float64})
    vector = deepcopy(vec_input)
    inv_matrix = deepcopy(input)
    rows, collumns = size(inv_matrix)
    if rows != collumns
        error("Matrix is not squared")
    end
    for a in (1:rows-1)
        #partial pivot
        max = argmax(abs.(inv_matrix[a:rows, a]))+a-1
        # changing the rows
        tmpvec = inv_matrix[a,:]
        inv_matrix[a,:] = inv_matrix[max,:] 
        inv_matrix[max,:] = tmpvec
        tmpval = vector[a]
        vector[a] = vector[max]
        vector[max] = tmpval
        # until here
        factor = inv_matrix[a,a] 
        for i in (1:rows)
            inv_matrix[a,i] = inv_matrix[a,i] / factor
        end
        vector[a] = vector[a] / factor
        for c in (a+1:rows)            
            subtraction_factor = inv_matrix[c,a]
            for i in (1:rows)
                inv_matrix[c,i] = inv_matrix[c,i]-subtraction_factor*inv_matrix[a,i] 
            end
            vector[c] = vector[c]-subtraction_factor*vector[a]
        end
    end
    return inv_matrix, vector
end

function backward_substitution(input::Matrix, vector::Vector) #passt so
    upper_matrix  = deepcopy(input)
    vector2 = deepcopy(vector)
    #ab hier die substitution
    rows, collumns = size(upper_matrix)
    out = Vector{Float64}(undef, rows)
    for x in LinRange(rows,1, rows)
        subtraction = 0
        for i in (x+1:rows)
            subtraction += upper_matrix[x,i]*out[i] 
        end
        out[x] = (vector2[x]-subtraction)/upper_matrix[x,x]
    end
    return out
end

function calculate_inverse(A::Matrix)
    rows, collumns = size(A)
    if rows != collumns
        error("Matix not squared")
    end
    A_inverse = Matrix{Float64}(undef, rows, collumns)
    for r in (1:rows)
        b = Float64[i==r for i in (1:rows)]
        test, x = GaußJordan_w_pivot(A, b)
        A_inverse[:,r] = backward_substitution(test, x)
    end
    return A_inverse
end


#TODO: Aufgabe
#TODO: pivot einrichten, sonst geht das nicht wegen null an a,a stelle
aufgabenmatrix = Float64[3 -1 1 4 1;1 0 1 1 1;-1 -1 -1 -1 -1;0 0 0 7 2;3 4 -3 4 1]
b0 = Float64[1,-1,1,-1,1]
b1 = Float64[0,1,2,3,4]


triag, y0 = GaußJordan_w_pivot(aufgabenmatrix, b0)
x0= backward_substitution(triag, y0)
aufgabenmatrix*x0
triag, y1 = GaußJordan_w_pivot(aufgabenmatrix, b1)
x1 = backward_substitution(triag, y1)
aufgabenmatrix*x1


#Warum geht die inverse nicht mehr???
a_inverse = calculate_inverse(aufgabenmatrix)
aufgabenmatrix*a_inverse



include("../Matrixfunctions.jl")
In = [1.1 1.7 1.3;1.3 1.9 2.3; 2.1 3.1 2.9]
test = calculate_inverse(In)
dx0 = test * [-0.004,-0.076,-0.088]
x0 = [-20.1,13.4,1.02]
b = [2.0,1.6,2.2]

test
x0+dx0
In*(x0+dx0)
In*x0
b-In*x0