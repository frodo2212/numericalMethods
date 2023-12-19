

function betrag(x)
    return sqrt(sum(map(xi->xi^2,x)))
end


function power(A, x0; error = 1e-10, maxstep = 1000000)
    x1 = x0
    lambda1 = 0
    for i in (1:maxstep)
        xn = A*x1
        x1 = xn/betrag(xn)
        lambda1 = x1'*(A*x1)
        if betrag(A*x1-lambda1*x1) < error
            break
        end
    end
    # if !isapprox(A*x1,lambda1*x1)
    #     error("Problem mit 1. Eigenvector")
    # end
    
    x2 = x0-(x0'*x1)*x1
    lambda2 = 0
    for i in (1:maxstep)
        xp2 = A*x2
        alpha = xp2'*x1
        x2 = (xp2-alpha*x1)/betrag(xp2-alpha*x1)
        lambda2 = x2'*(A*x2)
        if betrag(A*x2-lambda2*x2) < error
            break
        end
    end
    # if !isapprox(A*x2,lambda2*x2)
    #     error("Problem mit 2. Eigenvector")
    # end

    return x1,x2, lambda1,lambda2
end

A = [1 -2 2 ; -2 0 0 ; 2 0 2]
x0 = [1,1,1]

x1,x2, l1,l2 = power(A,x0,error=1e-13)
