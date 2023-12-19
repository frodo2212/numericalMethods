function myphia(x)
    return exp(-x)
end 
function myphib(x)
    return (1+x)/(1+exp(x))
end 
function myphic(x)
    return 1+x-x*exp(x)
end 
function myphiNR(x)
    # return x - (x*exp(x)-1)/(exp(x)*(1+x))
    return (x^2+exp(-x))/(1+x) 
end 

function testf(x)
    return x*exp(x)-1
end

function W71(x)
    return x-(x^3-5*x-1)/(3*x^2-5)
end
function W72(x)
    return x+(cos(x)-x)/(sin(x)+1)
end


function NRmethod(phi, x0, n)
    values = Vector{Float64}(undef, n)
    values[1] = phi(x0)
    for i in (2:n)
        values[i] = phi(values[i-1])
    end
    return values
end

testa = NRmethod(myphia, 0.5, 10)
testb = NRmethod(myphib, 0.5, 10)
testc = NRmethod(myphic, 0.5, 10)
testNR = NRmethod(myphiNR, 0.5, 10)


#für die neue woche
neu1 = NRmethod(W71,1,5)
neu2 = NRmethod(W72,1,5)
cos(0)*(cos(2)-2)
a=2
1*(a^3-5*a+1)




testf(testNR[40])

testall(n=1)

function testall(;x0=0.5, n=10)
    testa = NRmethod(myphia, x0,n)
    testb = NRmethod(myphib, x0,n)
    testc = NRmethod(myphic, x0,n)
    testNR = NRmethod(myphiNR, x0,n)
    println(string("function a results in:              ",testf.(testa[n])))
    println(string("function b results in:              ",testf.(testb[n])))
    println(string("function c results in:              ",testf.(testc[n])))
    println(string("Newton-Raphson-function results in: ",testf.(testNR[n])))
end
