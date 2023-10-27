function tridiag(input::Matrix, vector::Vector)
    rows, collumns = size(input)
    y = Vector{Float64}(undef, rows)
    x = Vector{Float64}(undef, rows)
    out = Matrix{Float64}(undef, rows, collumns)
    # alpha = Matrix{Float64}(undef, rows, collumns)
    # beta = Matrix{Float64}(undef, rows, collumns)
    y[1] = vector[1]
    out[1,1] = input[1,1]
        # beta[1,1]=input[1,1]
        # alpha[1,1] = 1
    for i in (2:rows)
        out[i-1,i]=input[i-1,i]
            # beta[i-1,i]=input[i-1,i]
        out[i,i-1] = input[i,i-1]/out[i-1,i-1]
            # alpha[i,i-1] = out[i,i-1]
        out[i,i] = input[i,i]-out[i,i-1]*out[i-1,i]
            # beta[i,i] = out[i,i]
            # alpha[i,i] = 1
        y[i] = vector[i]-out[i,i-1]*y[i-1]
    end
    x[rows] = y[rows]/out[rows,rows]
    for i in [3,2,1] 
        x[i] = (y[i]-(out[i,i+1]*x[i+1]))/out[i,i]
    end
    return out, x#, y# , alpha, beta
end

tridiag_matrix = Float64[2 -1 0 0;-1 2 -1 0; 0 -1 2 -1; 0 0 -1 2]
testvec = [-5,1,4,1]
test, vec= tridiag(tridiag_matrix, testvec)
tridiag_matrix*vec



#alles ab hier ist irrelevant, das ist nur de crout ohne pivoting 


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


testmatrix = Float64[3 -1 1 4 1;1 1 1 3 1;0 -1 -1 -1 -1;0 0 0 7 2;1 4 -3 4 1]
testvec = Float64[1,-1,1,-1,1]
test = crout(testmatrix)
y = forward_substitution(test, testvec)
x = backward_substitution(test, y)
testmatrix*x