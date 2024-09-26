using BenchmarkTools

function test1(x, N)
    tmp = zeros(2)
    for n in 1:N
        tmp[1] = n * rand()
        tmp[2] = n * rand()
        x[n, :] .= tmp
    end
end

function test2(x::Vector)
    
end

x = zeros(100, 2)
@btime test1(x, 10)