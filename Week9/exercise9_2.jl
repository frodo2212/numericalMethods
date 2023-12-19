matrix1 = Float64[-1 3 2 1;3 -3 1 5;2 1 2 -4;1 5 -4 1]
matrix2 = zeros(Float64, 100,100)
[matrix2[i,i] = 2 for i in (1:100)]
[matrix2[i-1,i] = -1 for i in (2:100)]
[matrix2[i,i-1] = -1 for i in (2:100)]
matrix2
#Formeln auf seite 109 im skript

using LinearAlgebra
eigvals(matrix1)
eigvals(matrix2)


m1, ev = transform(matrix1, output=true, showstuff=true)
test = sum(diag)
diag = [matrix1[i,i] for i in (1:4)]
a = 4
matrix1*vec(ev[:,a])
ev[:,a]*m1[a,a]
ev



ev[:,1]


transform(matrix2)
test = sum(diag)
diag = [matrix2[i,i] for i in (1:100)]
matrix2

function transform(umatrix::Matrix; maxsweeps=50, err=1e-5, verr=1e-5, showstuff::Bool=false, output::Bool=false)
    matrix = deepcopy(umatrix)
    dim = size(matrix)[1]
    ev = zeros(Float64, dim, dim)
    [ev[i,i] = 1 for i in (1:dim)]
    rot = 0
    sweeps = 0
    for sweep in (1:maxsweeps)
        sweeps = sweep
        sum = 0
        for i in (1:dim)
            for j in (i+1:dim)
                sum += abs(matrix[i,j])
            end
        end
        if sum < err
            println("zero offdiagonal elementes reached in sweep $sweep")
            break
        end
        for i in (1:dim)
            for j in (i+1:dim)
                if !isapprox(matrix[i,j],0)
                    rot += 1
                    jacobi1(dim, i, j, matrix, ev)
                end
            end
        end
    end
    if sweeps == maxsweeps
        println("maximum number of sweeps reached!")
    end
    println("found Eigenvalues and Eigenvectors after $rot Jacobi-rotations")
    everr = Int64[]
    for i in (1:dim)
        lambda = matrix[i,i]
        evi = vec(ev[:,i])
        diff = [(umatrix*evi)[j]-(lambda*evi)[j] for j in (1:dim)]
        if sum(abs.(diff)) < verr*dim
            showstuff && println("lambda[$i] = $lambda")
            showstuff && print("Eigenvector x[$i] = ")
            showstuff && println(evi)
        else
            showstuff && println("Eigenvector x[$i] is wrong")
            push!(everr, i)
        end
    end
    if everr == Int64[]
        println("everything sucessful")
    else
        print("Eigenvectors which did not work")
        println(everr)
    end
    output && return matrix, ev
end



function jacobi1(dim, p, q, matrix, ev)    
    pqalt = matrix[p,q]
    ppalt = matrix[p,p]
    qqalt = matrix[q,q]
    # TODO: irgendwie theta errechnen
    vartheta = (qqalt-ppalt)/(2*pqalt)
    t = sign(vartheta)/(abs(vartheta)+sqrt(vartheta^2+1))
    if abs(atan(t)) > pi/4
        t = 1/(2*vartheta)
        print("overflow")
    end
    if t == 0
        t = 1
    end
    c = 1/sqrt(t^2+1)
    s = t/sqrt(t^2+1)
    matrix[p,p] = c*c*ppalt+s*s*qqalt-2*s*c*pqalt
    matrix[q,q] = c*c*qqalt+s*s*ppalt+2*s*c*pqalt
    matrix[p,q] = (c*c-s*s)*pqalt+s*c*(ppalt-qqalt)
    matrix[q,p] = matrix[p,q]
    for r in (1:dim)
        if r != p && r != q
            rqalt = matrix[r,q]
            rpalt = matrix[r,p]
            matrix[r,q] = c*rqalt+s*rpalt
            matrix[q,r] = matrix[r,q]
            matrix[r,p] = c*rpalt-s*rqalt
            matrix[p,r] = matrix[r,p]
        end
    end
    #trafo von den EVs
    #s = -s
    for r in (1:dim)
        rqalt = ev[r,q]
        rpalt = ev[r,p]
        ev[r,q] = c*rqalt+s*rpalt
        ev[r,p] = c*rpalt-s*rqalt
    end
end

function jacobi(dim, p, q, matrix)
    #(app-aqq)/apq = (c^2+s^2)/(c*s)
    # TODO: irgendwie theta errechnen
    vartheta = (app-aqq)/(2*apq)
    t = sgn(vartheta)/(abs(vartheta)+sqrt(vartheta^2+1))
    # theta = arctan(t) #- brauch ich das überhaupt?
    # if irgendwas zu groß
    #     t = 1/(2*vartheta)
    # end
    c = 1/sqrt(t^2+1)
    s = t/sqrt(t^2+1)
    Jac = zeros(Float64, dim, dim)
    [Jac[i,i] = 1 for i in (1:dim)]
    Jac[p,q] = -s #aufpassen wo hier das minus hin muss...
    Jac[p,p] = c
    Jac[q,p] = s
    Jac[q,q] = c
    return Jac
end
