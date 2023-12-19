using CairoMakie

function test(q::Float64 ;params=(pi,1,1,0.5)) #a,m,k1,k2
    tmp = params[3]*sin(q*params[1]/2)^2+params[4]*sin(q*params[1])^2
    return 2/sqrt(params[2])*sqrt(tmp)
end

function plotIt(f, intervall; samples=500)
    fig = Figure()
    ax = Axis(fig[1,1])
    xtest = LinRange(intervall[1],intervall[2],samples)
    lines!(ax, xtest, f.(xtest, params=(pi,4,1,0.1)), label="k2=0.1")
    lines!(ax, xtest, f.(xtest, params=(pi,4,1,0.3)), label="k2=0.3")
    lines!(ax, xtest, f.(xtest, params=(pi,4,1,0.5)), label="k2=0.5")
    lines!(ax,[-1,-1],[0,1.1], linestyle=:dash, color=:darkgreen)
    lines!(ax,[1,1],[0,1.1], linestyle=:dash, color=:darkgreen)
    axislegend()
    return fig
end

xtest = LinRange(0,1,10)
test.(xtest, params=(pi,1,1,0))

fig = plotIt(test, (-2,2))
fig
save("EP-Bild.png",fig)


#
# f(x) = x^3-5*x+1
# fig = Figure()
# ax = Axis(fig[1,1])
# xtest = LinRange(0,2,100)
# lines!(ax,xtest,f.(xtest))
# fig


function plot3D(f;intx=(0,5),inty=(0,5),intz=(0,5),samples=50)
    fig = Figure()
    ax = Axis3(fig[1,1])
    xtest = LinRange(intx[1],intx[2],samples)
    ytest = LinRange(inty[1],inty[2],samples)
    ztest = LinRange(intz[1],intz[2],samples)
    z = [f(x,y) for x in xtest, y in ytest]
    #scatter!(xtest,ytest,z)
    #contour!(ax,xtest, ytest, f2)
    #contour3d!(xtest, ytest, z, linewidth=2, color=:red2)
    surface!(ax,xtest, ytest, z)
    return fig
end
f(x,y) = 2*sin(x)+x*cos(y)
f2(x,y) = x*y^3-y*x^3

plot3D(f)
plot3D(f2)
