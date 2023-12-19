include("../Matrixfunctions.jl")
# using LinearAlgebra

function fnkt(x1,x2,x3)
    return x1^2+x2-11, x1+x2^2-7, x1*x2-x3^3
end
function myfunction(x::Vector)
    y = fnkt(x[1],x[2],x[3])
    Jac = Matrix{Float64}(undef, 3,3)
    Jac[1,1] = 2*x[1]
    Jac[2,1] = 1
    Jac[3,1] = x[2]
    Jac[1,2] = 1
    Jac[2,2] = 2*x[2]
    Jac[3,2] = x[1]
    Jac[1,3] = 0
    Jac[2,3] = 0
    Jac[3,3] = -3*x[3]^2
    return y, Jac
end

function solve(F, Jac::Matrix)
    Fvec = [i for i in F]
    delta = calculate_inverse(Jac)*Fvec
    return delta
end

function f(x::Vector)
    y = fnkt(x[1],x[2],x[3])
    return 0.5*(y[1]^2+y[2]^2+y[3]^2)
end

function fprime(x, Jac,  deltax)
    return (Matrix(Jac')*[i for i in fnkt(x[1],x[2],x[3])])'*deltax
end

function main(x0::Vector; tolf = 1e-6, tolx = 1e-6, MAXSTEPS = 200, showinfo::Bool=false, out::Bool=false)
    println("main 1: using start-value: $x0")
    root=x0
    Step = 0
    for step in (1:MAXSTEPS)
        F, Jac = myfunction(root)
        errf = sqrt(F[1]^2+F[2]^2+F[3]^2)
        errf < tolf && break # exit loop if small step
        showinfo && print(string("step $step: errf = ",round(errf,digits=6)))
        deltax = solve(F,Jac)
        root -= deltax
        Step += 1
        errx = sqrt(deltax[1]^2+deltax[2]^2+deltax[3]^2)
        showinfo && println(string(", errx = ",round(errx,digits=6)))
        errx < tolx && break # exit loop if small step
    end
    println(string(Step, " steps used"))
    print("   x =  $root\n")
    solution = fnkt(root[1],root[2],root[3])
    println("f(x) = $solution\n")
    out && return root, Step
end

function main2(x0::Vector; tolf = 1e-6, tolx = 1e-6, MAXSTEPS = 200, showinfo::Bool=false, out::Bool=false)
    println("main 2: using start-value: $x0")
    root=x0
    Step = 0
    for step in (1:MAXSTEPS)
        F, Jac = myfunction(root)
        errf = sqrt(F[1]^2+F[2]^2+F[3]^2)
        errf < tolf && break # exit loop if small step
        showinfo && print(string("step $step: errf = ",round(errf,digits=6)))
        deltax = solve(F,Jac)
        #TODO: Hier ist das neue Zeug drin
        if f(root-deltax) < f(root)
            root-=deltax
        else
            for i in (1:20)
                if f(root-deltax*(20-i)/20) < f(root)
                    root -= deltax*(20-i)/20
                    println(string("stepsize reduced to ", (20-i)/20, " in step ", step))
                    break
                end
                if i == 20
                    # println(string("iterationen scheiße",step))
                    root -= deltax*0.01
                end
            end
        end
        Step += 1
        errx = sqrt(deltax[1]^2+deltax[2]^2+deltax[3]^2)
        showinfo && println(string(", errx = ",round(errx,digits=6)))
        errx < tolx && break # exit loop if small step
    end
    println(string(Step, " steps used"))
    print("   x =  $root\n")
    solution = fnkt(root[1],root[2],root[3])
    println("f(x) = $solution\n")
    out && return root, Step
end

function main3(x0::Vector; tolf = 1e-6, tolx = 1e-6, MAXSTEPS = 200, showinfo::Bool=false, out::Bool=false)
    println("main 3: using start-value: $x0")
    root=x0
    Step = 0
    for step in (1:MAXSTEPS)
        F, Jac = myfunction(root)
        errf = sqrt(F[1]^2+F[2]^2+F[3]^2)
        errf < tolf && break # exit loop if small step
        showinfo && print(string("step $step: errf = ",round(errf,digits=6)))
        deltax = solve(F,Jac)
        #TODO: Hier ist das neue Zeug drin
        if f(root-deltax) < f(root)
            root-=deltax
        else
            lambda_alt = 1
            for i in (1:20)
                g0p = fprime(root,Jac,deltax)
                g1 = f(root+deltax*lambda_alt)
                g0 = f(root)
                lambda = -(g0p)/(2*(g1-g0-g0p))
                if lambda < 0.1*lambda_alt #Muss das hier 0.1*lambda_alt sein oder nur 0.1?
                    lambda = 0.1*lambda_alt
                end
                if f(root-deltax*lambda) < f(root)
                    root -= deltax*lambda
                    println(string("stepsize reduced to ", lambda, " in step ", step," with ",i," iterations"))
                    break
                end
                lambda_alt = lambda
                if i == 20
                    println(string("iterationen scheiße ",step," root at ", root))
                    root -= deltax*0.01
                end
            end
        end
        Step += 1
        errx = sqrt(deltax[1]^2+deltax[2]^2+deltax[3]^2)
        showinfo && println(string(", errx = ",round(errx,digits=6)))
        errx < tolx && break # exit loop if small step
    end
    println(string(Step, " steps used"))
    print("   x =  $root\n")
    solution = fnkt(root[1],root[2],root[3])
    println("f(x) = $solution\n")
    out && return root, Step
end


main([1,1,1], showinfo=true, out=true)
main([-1,-1,-1], showinfo=true, out=true)
main([0.1,0.1,0.1], showinfo=true, out=true)
main([-0.1,-0.1,-0.1], showinfo=true, out=true)
# main([0,0,0])



main2([1,1,1], tolf=1e-12)
main3([1,1,1])

function compare(x;tolf = 1e-6, tolx = 1e-6)
    r1, step1 = main(x, tolf=tolf, tolx=tolx, out=true)
    r2, step2 = main2(x, tolf=tolf, tolx=tolx, out=true)
    r3, step3 = main3(x, tolf=tolf, tolx=tolx, out=true)
    println(string("method 1 took ", step1," steps and found root: ",r1))
    println(string("method 2 took ", step2," steps and found root: ",r2))
    println(string("method 3 took ", step3," steps and found root: ",r3))
end
compare([1,1,1], tolx=1e-12, tolf=1e-12)


compare([1,1,1])
compare([-1,-1,-1])
compare([0.1,0.1,0.1])
compare([-0.001,-0.001,-0.001])