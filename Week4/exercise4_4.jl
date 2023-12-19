using CairoMakie
include("../Matrixfunctions.jl")

function P_params(x_in, y_in, m, i)
    if m == 0
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
    return p_new
end

function generateMatrix(n, x, y)
    #Hier wird im moment implizit 0 als randbedingungen angenommen
    matrix = zeros(Float64, n-2, n-2)
    [matrix[j,j] = (x[j+2]-x[j])/3 for j in (1:n-2)]
    [matrix[j,j-1] = (x[j]-x[j-1])/6 for j in (2:n-2)]
    [matrix[j-1,j] = (x[j+1]-x[j])/6 for j in (2:n-2)]
    b = [(y[j+1]-y[j])/(x[j+1]-x[j])-(y[j]-y[j-1])/(x[j]-x[j-1]) for j in (2:n-1)]
    return matrix, b
end

function cubic_spline(xdata::Vector{T}, ydata::Vector{T}; showinfo::Bool=false) where T<:Real
    len = length(xdata)
    if len != length(ydata)
        error("sizes of the two input vectors dont match")
    end
    matrix, b = generateMatrix(len, xdata, ydata)
    yss = Float64[0]
    append!(yss, tridiag_x_calculation(matrix, b))
    push!(yss, 0)
    showinfo && @show yss

    A(x, j) = (xdata[j+1]-x)/(xdata[j+1]-xdata[j])
    B(x, j) = (x-xdata[j])/(xdata[j+1]-xdata[j])
    ylin(x, j) = A(x, j)*ydata[j]+B(x, j)*ydata[j+1]
    C(x, j) = (1/6)*((A(x, j))^3-A(x, j))*((xdata[j+1]-xdata[j])^2)
    D(x, j) = (1/6)*((B(x, j))^3-B(x, j))*((xdata[j+1]-xdata[j])^2)
    pfn(x, j) = C(x, j)*yss[j]+D(x, j)*yss[j+1]
    ynew(x,j) = ylin(x,j)+pfn(x,j)
    #ys(x,j) = (ydata[j+1]-ydata[j])/(xdata[j+1]-xdata[j])-(3*A(x,j)^2-1)/(6)*(xdata[j+1]-xdata[j])*yss[j]+(3*B(x,j)^2-1)/(6)*(xdata[j+1]-xdata[j])*yss[j+1]

    for j in (1:len-1) 
        if !isapprox(pfn(xdata[j],j),0) || !isapprox(pfn(xdata[j+1],j),0)
            error("polynom erfüllt die Randbedingungen nicht")
        end
        if !isapprox(ylin(xdata[j+1],j),ydata[j+1]) || !isapprox(ylin(xdata[j],j),ydata[j])
            error("ylin interpolation hat nicht funktioniert")
        end
    end
    function yfinished(x)
        y=0
        if x < xdata[1] || x > xdata[len]
            error("value out of bound")
        end
        for j in (1:len-1)
            if xdata[j] <= x && xdata[j+1] >= x
                y = ynew(x,j)
            end
        end
        return y
    end
    return yfinished
end

function plot_cubic_spline(xdata, ydata; addpoints::Bool=true, addlinear::Bool=false)
    yinterpol = cubic_spline(xdata, ydata)
    fig = Figure()
    ax = Axis(fig[1,1])
    addpoints && scatter!(ax, xdata, ydata)
    steps = LinRange(xdata[1],xdata[length(xdata)],120)
    lines!(ax, steps, yinterpol.(steps), color = :darkred)
    if addlinear 
        for i in (1:length(xdata)-1)
            lines!(ax, [xdata[i],xdata[i+1]], [ydata[i],ydata[i+1]], linestyle = :dash, color = :green)
        end
    end
    return fig
end


xdata = [1.1, 1.9, 3.2, 4.1, 4.9, 6.3]
ydata = [6.0, 5.1, 3.9, 2.2, 4.0, 6.2]
xdata2 = [0.9, 1.8, 3.1, 3.9, 4.8, 6.1, 7]
ydata2 = [1.0, 5.5, 6.0, 1.4, 0.9, 4.5, 4.5]


fig = plot_cubic_spline(xdata, ydata, addlinear=true)
fn = cubic_spline(xdata, ydata)
fn(4)
fn(1.9)





#Aufgabe 2
xdataA2 = Float64[0,0.5,1,2]
ydataA2 = Float64[1,0.8,0.5,0.2]

fn = cubic_spline(xdataA2, ydataA2, showinfo=true)
fig = plot_cubic_spline(xdataA2, ydataA2)