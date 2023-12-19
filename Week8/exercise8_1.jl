# using CairoMakie

function mypoly2(x)
    return (((((((x-3.91823)*x-7.07373)*x+39.9139)*x-5.12949)*x-108.519)*x+72.822)*x+72.5235)*x-61.688
end

function polynome(x, a)
    len = length(a)
    y = 0
    for i in (1:len-1)
        y = (y+a[i])*x
    end
    return y+a[len]
end


function df(a0::Vector{Float64})
    len = length(a0)
    a1 = Vector{Float64}(undef, len-1)
    for i in (1:len-1)
        a1[i] = a0[i]*(len-i)
    end
    return a1
end


a = [1.0,-3.91823,-7.07373,39.9139,-5.12949,-108.519,72.822,72.5235,-61.688]

function mayhly(a, xguess0::Vector{Float64}; xerr=1e-14, ferr=1e-12, Maxstep=300)
    xguess = deepcopy(xguess0)
    aprime = df(a)
    for j in (1:length(xguess))
        for i in (1:Maxstep)
            y = polynome(xguess[j],a)
            fs = polynome(xguess[j], aprime)
            nenner = [1/(xguess[j]-xguess[q]) for q in (1:j-1)]
            pnj = y*sum(nenner)
            dx = y/(fs-pnj)
            xguess[j] = xguess[j]-dx
            if abs(dx) < xerr
                break
            end 
        end
    end
    if sum(abs.(mypoly2.(xguess))) > ferr*length(xguess)
        error("some value is non zero")
    end
    print("polished x values: ")
    println(round.(xguess,digits=3))
    print("respective function values: ")
    println(round.(mypoly2.(xguess),digits=3))
    return xguess
end



x0 = [0.9,-1.1,1.4,1.4,-2.0,-2.0,2.2,2.2]
test1 = mayhly(a,x0)
test1
mypoly2.(test1)
mypoly2.(x0)


function plot(a, int, ylim=(0,0))
    fig = Figure()
    ax = Axis(fig[1,1])
    xsteps = LinRange(int[1],int[2],1000000)
    yvalues = [polynome(x,a) for x in xsteps]
    yvalues = [mypoly2(x) for x in xsteps]
    if ylim != (0,0)
        ylims!(ylim[1],ylim[2])
    end
    lines!(ax, xsteps, yvalues)
    return fig
end
fig1 = plot(a,(-10,10))
fig2 = plot(a, (-4,6))
fig3 = plot(a,(-2.6,3.2))
fig4 = plot(a,(0,1.5))

save("Polynom.png",fig3)