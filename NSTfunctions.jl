include("Matrixfunctions.jl")

function bisection(f, x1, x2, xacc)
    if myfunction(x2)*myfunction(x1)>0 
        println("root not bracketed!")
    end
    for i in (1:300)
        middle = (x1+x2)/2
        if f(middle) == 0 || x2-x1 < xacc            
            println(string("bisection    root found after ",i,"steps at: ",middle))
            return middle
        elseif f(x1)*f(middle) < 0
            x2 = middle
        elseif f(x2)*f(middle) < 0
            x1 = middle
        else
            println("Problemchen")
            return 0
        end
    end
    println("maximum amount of steps reached")
    return (x1,x2)
end

function secant(fkt, x1, x2, xacc)
    xneg = 0.0
    root = 0.0
    f1=fkt(x1)
    f2=fkt(x2)
    if f1*f2>0 
        println("root not bracketed!")
    end
    if f1 < 0 
        xl=x1
        root=x2
    else
        xl=x2
        root=x1
    end
    fneg = fkt(xneg)
    fpos = fkt(root)
    for i in (1:40)
        dx = (root-xneg)*(fpos)/(fneg-fpos)
        xneg = root
        fneg = fpos
        root += dx
        fpos = fkt(root)
        if abs(dx) < xacc || fpos == 0            
            println(string("secant       root found after ",i,"steps at: ",root))
            return root
        end
    end
    println("max steps reached")
    return root
end

function regulafalsi(f, x1, x2, xacc)
    fa = f(x2)
    fb = f(x1)
    if fa*fb>0 
        println("root not bracketed!")
    end
    a = x2
    b = x1
    old = (b-a)/2
    for i in (1:100)
        c = (a*fb-b*fa)/(fb-fa)
        fnew = f(c)
        if fnew*fa > 0
            old = a
            a = c
            fa = fnew
        else
            old = b
            b = c
            fb = fnew
        end
        if abs(old-a) < xacc || abs(old-b) < xacc || fnew == 0
            println(string("regulafalsi  root found after ",i,"steps at: ",a))
            return a
        end
    end
    println("maximum steps achieved")
    return b, a
end

function NRmethod(phi, x0, n)
    values = Vector{Float64}(undef, n)
    values[1] = phi(x0)
    for i in (2:n)
        values[i] = phi(values[i-1])
    end
    return values
end


