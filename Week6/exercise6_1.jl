X10=[-0.9739065285171717,-0.8650633666889845,-0.6794095682990244,-0.4333953941292472,-0.1488743389816312,0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717]
w10=[0.0666713443086881,0.1494513491505806,0.2190863625159821,0.2692667193099963,0.2955242247147529,0.2955242247147529,0.2692667193099963,0.2190863625159821,0.1494513491505806,0.0666713443086881]
X2=[-sqrt(1/3),sqrt(1/3)]
w2=[1,1]
X3=[-sqrt(3/5),0,sqrt(3/5)]
w3=[10/18,8/9,10/18]
X4=[-sqrt((15+2*sqrt(30))/(35)),-sqrt((15-2*sqrt(30))/(35)),sqrt((15-2*sqrt(30))/(35)),sqrt((15+2*sqrt(30))/(35))]
w4=[0.717161239,0.26850642,0.26850642,0.717161239]

include("../Week5/exercise5_3.jl")

function myfunction(x)
    return sqrt(1+x^2)
end

function integrate(fkt, x1, x2, x, w)
    # intervall ändern
    f(x) = fkt(x*(x2-x1)/2+(x1+x2)/2)
    integral = 0
    len = length(x)
    if len != length(w)
        error("dimensionen stimmen nicht")
    end
    for i in (1:len) 
        integral += w[i]*f(x[i])
    end
    return (x2-x1)/2*integral
end

x1 = 0
x2 = 1
xacc = 1e-18

value2 = integrate(myfunction, x1, x2, X2, w2)
value3 = integrate(myfunction, x1, x2, X3, w3)
value4 = integrate(myfunction, x1, x2, X4, w4)
value10 = integrate(myfunction, x1, x2, X10, w10)
test1 = quadrature_trapez(myfunction, x1,x2,xacc)
test2 = quadrature_simpson(myfunction, x1,x2,xacc,maxstep=28)
test3 = quadrature_romberg(myfunction, x1,x2,xacc)

testfehler = abs(test3-test2)
fehler = abs(value10-test2)