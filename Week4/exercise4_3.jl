using CairoMakie


function sortArgs(xdata, ydata)
    len = length(xdata)
    if len != length(ydata)
        error("dimensions of x and y data dont match")
    end
    sortedx = deepcopy(xdata)
    sortedy = deepcopy(ydata)
    for z in (1:len)
        for i in (1:len-z)
            if sortedx[i] > sortedx[i+1]
                tmp = sortedx[i]
                sortedx[i] = sortedx[i+1]
                sortedx[i+1] = tmp

                tmp = sortedy[i]
                sortedy[i] = sortedy[i+1]
                sortedy[i+1] = tmp
            end
        end
    end
    return sortedx, sortedy
end

function locateX(x::Real, xdata::Vector)
    len = length(xdata)
    jlow = 1
    for i in (1:len-1)
        if xdata[i] <= x && xdata[i+1] >= x
            jlow = i
            break
        end
    end
    return jlow
end

function P_params(x_in, y_in, m::Int64, i::Int64; returnfunction::Bool=false)
    if m == 0
        returnfunction && (return f(x)=y_in[i])
        return [y_in[i]]
    end
    p1 = P_params(x_in,y_in,m-1,i)
    p2 = P_params(x_in,y_in,m-1,i+1)
    p_new = Float64[]
    push!(p_new, (p1[1] - p2[1])/(x_in[i]-x_in[i+m]))
    for j in (1:m-1)
        push!(p_new, (p1[j+1]-p2[j+1]+x_in[i]*p2[j]-x_in[i+m]*p1[j])/(x_in[i]-x_in[i+m]))
    end
    push!(p_new, (x_in[i]*p2[m]-x_in[i+m]*p1[m])/(x_in[i]-x_in[i+m]))
    if returnfunction
        function generate_y(xvalue::Real)
            order = length(p_new)
            yinter = 0
            for m in (1:order)
                yinter += p_new[m]*xvalue^(order-m)
            end
            return yinter
        end
        return generate_y 
    end
    return p_new
end



function plotIt(xdata, ydata)
    order = length(xdata)-1
    fig = Figure()
    ax = Axis(fig[1,1])
    scatter!(ax, xdata, ydata)
    yinterpol = P_params(xdata,ydata,order,1, returnfunction=true)
    xsteps = collect(LinRange(xdata[1],xdata[order+1],100))
    lines!(ax, xsteps, yinterpol.(xsteps))
    return fig
end

function plotIt(xdata, ydata, x, order; limitOrder::Bool=false)
    order -= 1 #macht das sinn so?
    xnew, ynew = sortArgs(xdata, ydata)
    if x < xnew[1] || x > xnew[length(xnew)]
        error("x not in the given xdata range")
    end
    len = length(xnew)
    fig = Figure()
    ax = Axis(fig[1,1])
    scatter!(ax, xnew, ynew)
    jlow = locateX(x, xnew)
    startpos = maximum([1,jlow-floor(Int64, order/2)])
    startpos + order > len && (startpos = maximum([1,len-order]))
    startpos + order > len && (order = len-startpos)
    if limitOrder
        if startpos==1
            order = (jlow-startpos)*2
        end
        if startpos+order == length(xnew)
            int = jlow+order-1
            order = int*2
            startpos = jlow-int
        end
    end
    scatter!(ax, xdata[startpos:startpos+order], ydata[startpos:startpos+order])    
    yinterpol = P_params(xdata,ydata,order,startpos, returnfunction=true)
    scatter!(ax,x,yinterpol(x))
    xsteps = collect(LinRange(xdata[startpos],xdata[startpos+order],75))
    lines!(ax, xsteps, yinterpol.(xsteps))
    return fig
end

function interpolate(xdata, ydata, x, order; showdata::Bool=true)
    order = order - 1
    xnew, ynew = sortArgs(xdata, ydata)
    len = length(xnew)
    #if x < xnew[1] || x > xnew[length(xnew)]
    #    error("x not in the given xdata range")
    #end
    jlow = locateX(x, xnew)
    startpos = maximum([1,jlow-floor(Int64, order/2)])
    (startpos + order) > len && (startpos = maximum([1,len-order]))
    (startpos + order) > len && (order = len-startpos)
    showdata && println(string("x lies in between ",jlow," and ",jlow+1))
    showdata && println(string("lowest index: ", startpos, ", highest index: ", startpos+order))
    yinterpol = P_params(xnew, ynew, order, startpos, returnfunction=true)
    y = yinterpol(x)
    return y
end

function plotall(xdata, ydata, order)
    xnew,ynew = sortArgs(xdata,ydata)
    xsteps = LinRange(xnew[1],xnew[length(xdata)],300)
    interpoly = [interpolate(xnew,ynew,i,order, showdata=false) for i in xsteps]
    fig = Figure()
    ax = Axis(fig[1,1])
    scatter!(ax, xnew, ynew)
    lines!(ax,xsteps,interpoly)
    return fig
end


xdata = [1.1, 1.9, 3.2, 4.1, 4.9, 6.3]
ydata = [6.0, 5.1, 3.9, 2.2, 4.0, 6.2]
xdata2 = [0.9, 1.8, 3.1, 3.9, 4.8, 6.1]
ydata2 = [1.0, 5.5, 6.0, 1.4, 0.9, 4.5]

interpolate(xdata, ydata, 4, 3)
fig = plotIt(xdata, ydata, 2, 4)
fig = plotIt(xdata2, ydata2)


fig = plotall(xdata,ydata,10)


#A2 data
fig = plotIt(Float64[0,0.5,1,2],Float64[1,0.8,0.5,0.2],4,1)

