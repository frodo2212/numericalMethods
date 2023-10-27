function crout(input::Matrix) #without any pivoting
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

testmatrix = Float64[3 -1 1 4 1;1 1 1 3 1;0 -1 -1 -1 -1;0 0 0 7 2;1 4 -3 4 1]
testvec = Float64[1,-1,1,-1,1]
test = crout(testmatrix)

huhu = forward_substitution(test, testvec)
ergebnis = backward_substitution(test, huhu)

testmatrix*ergebnis

function forward_substitution(input::Matrix, vector::Vector)
    upper_matrix  = deepcopy(input)
    vector2 = deepcopy(vector)
    #ab hier die substitution
    rows, collumns = size(upper_matrix)
    out = Vector{Float64}(undef, rows)
    for x in (1:rows)
        subtraction = 0
        for i in (1:x-1) #TODO: die Summe umschreiben
            subtraction += upper_matrix[x,i]*out[i] 
        end
        out[x] = vector2[x]-subtraction
    end
    return out
end

function backward_substitution(input::Matrix, vector::Vector)
    upper_matrix  = deepcopy(input)
    vector2 = deepcopy(vector)
    #ab hier die substitution
    rows, collumns = size(upper_matrix)
    out = Vector{Float64}(undef, rows)
    for x in [5,4,3,2,1] #TODO: über nen linspace in negativer reihenfolge ersetzen
        subtraction = 0
        for i in (x+1:rows)
            subtraction += upper_matrix[x,i]*out[i] 
        end
        out[x] = (vector2[x]-subtraction)/upper_matrix[x,x]
    end
    return out
end