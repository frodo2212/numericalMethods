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
    xdata[len] < x && (jlow = len)
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

function interpolate(xdata, ydata, x, order; showdata::Bool=false)
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

