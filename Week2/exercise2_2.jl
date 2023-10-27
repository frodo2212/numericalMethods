
aufgabenmatrix = Float64[3 -1 1 4 1;1 0 1 1 1;-1 -1 -1 -1 -1;0 0 0 7 2;3 4 -3 4 1]


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
        #new stuff
        #partial pivot
        max = argmax(abs.(inv_matrix[a:rows, a]))+a-1
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

function backward_substitution(input::Matrix, vector::Vector; do_GaußJordan::Bool=false, with_pivot::Bool=true) #passt so
    if do_GaußJordan && with_pivot
        upper_matrix, vector2 = GaußJordan_w_pivot(input, vector)
    elseif do_GaußJordan && !with_pivot
        upper_matrix, vector2 = GaußJordan_wo_pivot(input, vector)
    else
        upper_matrix  = deepcopy(input)
        vector2 = deepcopy(vector)
    end
    #ab hier die substitution
    rows, collumns = size(upper_matrix)
    out = Vector{Float64}(undef, rows)
    for x in [5,4,3,2,1]
        subtraction = 0
        for i in (x+1:rows)
            subtraction += upper_matrix[x,i]*out[i] 
        end
        out[x] = (vector2[x]-subtraction)/upper_matrix[x,x]
    end
    return out
end

testmatrix = Float64[3 -1 1 4 1;1 1 1 3 1;0 -1 -1 -1 -1;0 0 0 7 2;1 4 -3 4 1]
huhu, test_vec = GaußJordan_wo_pivot(testmatrix, testvec)
huhu
huhu_out = backward_substitution(huhu, test_vec)

huhu2, test_vec2 = GaußJordan_w_pivot(testmatrix, testvec)
huhu2
huhu2_out = backward_substitution(huhu2, test_vec2)

testmatrix*huhu_out

#TODO: Aufgabe
#TODO: pivot einrichten, sonst geht das nicht wegen null an a,a stelle
aufgabenmatrix = Float64[3 -1 1 4 1;1 0 1 1 1;-1 -1 -1 -1 -1;0 0 0 7 2;3 4 -3 4 1]
b0 = Float64[1,-1,1,-1,1]
b1 = Float64[0,1,2,3,4]

x0= backward_substitution(aufgabenmatrix,b0,do_GaußJordan=true)
aufgabenmatrix*x0
x1 = backward_substitution(aufgabenmatrix,b1,do_GaußJordan=true)
aufgabenmatrix*x1
aufgabenmatrix

a_in = calculate_inverse(aufgabenmatrix)
aufgabenmatrix*a_in
#Aufgabe ende

argmax(hihi[1:5, 3])

testvec = Float64[1,-1,1,-1,1]
hihi = Float64[1 2 3 4 5; 0 3 4 7 5; 0 0 1 1 2; 0 0 0 4 3;0 0 0 0 1]
GaußJordan_wo_pivot(hihi, testvec)
vec_hihi = backward_substitution(hihi, testvec)
manuell_vec_hihi = Float64[-2/3,1/3,0,-1,1]
hihi*vec_hihi




function calculate_inverse(A::Matrix)
    rows, collumns = size(A)
    if rows != collumns
        error("Matix not squared")
    end
    A_inverse = Matrix{Float64}(undef, rows, collumns)
    for r in (1:rows)
        b = Float64[i==r for i in (1:rows)]
        A_inverse[:,r] = backward_substitution(A, b, do_GaußJordan=true)
    end
    return A_inverse
end

# tests
# hihi = Float64[1 2 3 4 5; 0 3 4 7 5; 0 0 1 1 2; 0 0 0 4 3;0 0 0 0 1]
# testvec = Float64[1,-1,1,-1,1]
# manuell_vec_hihi = Float64[-2/3,1/3,0,-1,1] # zur manuellen überprüfung was raus kommen müsste
  manuell_vec_hihi == backward_substitution(hihi, testvec)
  hihi*backward_substitution(hihi, testvec) == testvec
