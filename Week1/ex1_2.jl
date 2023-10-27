function integral(n::Integer)
    if n==0
        return 1-cos(1)
    else
        return 1-n*(n-1)*integral(n-2)
    end
end

function integral(n::Integer, decimals)
    if n==0
        return 1-round(cos(1), digits=decimals)
    else
        return 1-n*(n-1)*integral(n-2, decimals)
    end
end

@show integral(8)
@show integral(8, 5)
@show integral(8, 3)
