function tridiag(input::Matrix, vector::Vector)
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
    for i in [3,2,1] 
        x[i] = (y[i]-(out[i,i+1]*x[i+1]))/out[i,i]
    end
    return out, x
end

tridiag_matrix = Float64[2 -1 0 0;-1 2 -1 0; 0 -1 2 -1; 0 0 -1 2]
testvec = [-5,1,4,1]
test, vec= tridiag(tridiag_matrix, testvec)
tridiag_matrix*vec
vec
