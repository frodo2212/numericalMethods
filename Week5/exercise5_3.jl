include("../Week4/exercise4_3V2.jl")
function myfunction(x)
    y = x^4*log(x+sqrt(1+x^2))
    return y
end


function calculate_Innew(myfunction, a, b, n)
    intervalls = 2^(n-1)
    deltaN = (b-a)/(intervalls)
    sum = 0
    for k in (0:intervalls-1)
        sum += myfunction(a+(k+0.5)*deltaN)
    end
    return sum*deltaN
end

function quadrature_trapez(myfunction, a,b,eps; maxstep=25)
    sum = Float64[]
    push!(sum,(b-a)/2*(myfunction(a)+myfunction(b)))
    # push!(sum,0.5*(calculate_Innew(myfunction,a,b,1)+sum[1]))
    for i in (2:maxstep)
        push!(sum,0.5*(calculate_Innew(myfunction,a,b,i-1)+sum[i-1]))
        if i > 5
            if sum[i] == 0.0 && sum[i-1] == 0.0 || abs(sum[i]-sum[i-1]) < eps*abs(sum[i-1]) 
                println(string("quadrature_trapez has stopped after",i,"steps"))
                return sum[i]
            end
        end
    end
    println("Maximum number of steps reached in quadrature_trapez_stop!");
    return sum[maxstep]
end


function quadrature_simpson(myfunction,a,b,eps;maxstep=25)
    In = Float64[]
    Iimproved = Float64[]
    push!(In,(b-a)/2*(myfunction(a)+myfunction(b)))
    for i in (2:maxstep)
        push!(In,0.5*(calculate_Innew(myfunction,a,b,i-1)+In[i-1]))
        push!(Iimproved,4/3*In[i]-1/3*In[i-1])
        if i > 5
            if Iimproved[i-1] == 0.0 && Iimproved[i-2] == 0.0 || abs(Iimproved[i-1]-Iimproved[i-2]) < eps*abs(Iimproved[i-2]) 
                println(string("quadrature_trapez has stopped after",i,"steps"))
                return Iimproved[i-1]
            end
        end
    end
    println("Maximum number of steps reached in quadrature_simpson_stop!");
    return Iimproved[maxstep-1]
end

function quadrature_romberg(myfunction,a,b,eps;maxstep=23)
    In = Float64[]
    hs = Float64[]
    rom_values = Float64[]
    push!(In,(b-a)/2*(myfunction(a)+myfunction(b)))
    push!(hs, b-a)
    for i in (2:maxstep)
        push!(In,0.5*(calculate_Innew(myfunction,a,b,i-1)+In[i-1]))
        push!(hs, (b-a)^2/(2^(2*i-2)))
        push!(rom_values,interpolate(hs,In,0,3))
        if i > 5
            if rom_values[i-1] == 0.0 && rom_values[i-2] == 0.0 || abs(rom_values[i-1]-rom_values[i-2]) < eps*abs(rom_values[i-2]) 
                println(string("quadrature_romberg has stopped after",i,"steps"))
                return rom_values[i-1]
            end
        end
    end
    println("Maximum number of steps reached in quadrature_trapez_stop!");
    return rom_values
end

a = 0
b = 2
eps = 1e-10


test1 = quadrature_trapez(myfunction,a,b,eps)
test2 = quadrature_simpson(myfunction,a,b,eps,maxstep=25)
test3 = quadrature_romberg(myfunction,a,b,eps)
