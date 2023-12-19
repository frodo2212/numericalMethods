function my_function(x)
    y = x^4*log(x+sqrt(1+x^2))
    return y
end

function circle(x)
    r = 1
    if x <= r
        return sqrt(r^2-x^2)
    end
    return 0
end

function f3D(x,y,z)
    return (z^2+(sqrt(x^2+y^2)-3)^2) <= 1
end

function f3DV2(x,y,z)
    if x >= 1 && x <= 4 && y >= -3 && y <= 4 && z >= -1 && z <= 1
        if (z^2+(sqrt(x^2+y^2)-3)^2) <= 1
            return 1
        end
    end
    return 0
end



function generatePointsnD(n::Integer, dims::Integer, ranges::Vector{Tuple{Float64,Float64}})
    if dims != length(ranges)
       error("input ranges dont match the dimension")
    end
    points = Vector{Tuple}(undef,n)
    for i in (1:n)
        random = rand(dims)
        rescaled = [ranges[i][1]+random[i]*(ranges[i][2]-ranges[i][1]) for i in (1:dims)]
        points[i] = Tuple(Float64(x) for x in rescaled)
    end
    return points
end


function monte_carlo1DV2(myfunction, n, intervall=(0,1))
    points = generatePointsnD(n,1,Tuple{Float64,Float64}[intervall])
    integral = 0
    int2 = 0
    for point in points
        integral += myfunction(point[1])
        int2 += myfunction(point[1])^2
    end
    Volume = (intervall[2]-intervall[1])
    error = Volume*sqrt((int2/n-(integral/n)^2)/n)
    return Volume*integral/n, error
end

function monte_carlo1D(myfunction, n, intervall=(0,1), range=(0,1))
    points = generatePointsnD(n,2,Tuple{Float64,Float64}[intervall,range])
    below = 0
    for point in points
        if point[2] <= myfunction(point[1])
            below += 1
        end
    end
    Volume = (intervall[2]-intervall[1])*(range[2]-range[1])
    return Volume*(below / n)
end

function mc3D_Body_small(densityfunction, n, ranges)
    points = generatePointsnD(n, 3, ranges)
    xges = 0
    yges = 0
    zges = 0
    below = 0
    for point in points
        value = densityfunction(point[1],point[2],point[3])
        if value != 0
            below+=1
            xges += point[1]#*value
            yges += point[2]#*value
            zges += point[3]#*value
        end
    end
    Vi = (ranges[1][2]-ranges[1][1])*(ranges[2][2]-ranges[2][1])*(ranges[3][2]-ranges[3][1])
    Volume = Vi*(below / n)
    return Volume, (xges/below,yges/below,zges/below)
end

function mc3D_Body(densityfunction, n, ranges)
    points = generatePointsnD(n, 3, ranges)
    xges = 0
    yges = 0
    zges = 0
    below = 0
    av = 0
    av2 = 0
    xg2 = 0
    yg2 = 0
    zg2 = 0
    for point in points
        value = densityfunction(point[1],point[2],point[3])
        if value!=0
            below+=1
            av += value
            av2 += value^2
            xges += point[1]*value
            yges += point[2]*value
            zges += point[3]*value
            xg2 += (point[1]*value)^2
            yg2 += (point[2]*value)^2
            zg2 += (point[3]*value)^2
        end
    end
    Vi = (ranges[1][2]-ranges[1][1])*(ranges[2][2]-ranges[2][1])*(ranges[3][2]-ranges[3][1])
    # Volume = Vi*(below / n)
    Volume = Vi*(av/n) #eigentlich wäre das hier mass
    error = Vi*sqrt(((av2/n)-(av/n)^2)/n)
    vec_error = [Vi*sqrt(((xg2/n)-(xges/n)^2)/n),Vi*sqrt(((yg2/n)-(yges/n)^2)/n),Vi*sqrt(((zg2/n)-(zges/n)^2)/n)]
    return Volume, (xges/below,yges/below,zges/below), error, vec_error
end
        
i1 = monte_carlo1D(circle, 100000, (0, 1))*4
i2 = monte_carlo1DV2(circle, 100000, (0.0, 1.0))
i2[1]*4

# i1 = monte_carlo1DV2(my_function, 100000, (0, 0.5))


ranges = Tuple{Float64,Float64}[(1,4),(-3,4),(-1,1)]
i3 = mc3D_Body(f3DV2, 10000000, ranges)


println(string("Volume:",i3[1],"+-",i3[3]))
for i in (1:3)
    dir = ["x","y","z"]
    println(string("Center of mass ",dir[i],": ",i3[2][i],"+-",i3[4][i]))
end

i4 = mc3D_Body_small(f3D, 1000000, ranges)