using CSV
using DataFrames
using Infiltrator
using JLD

Base.@kwdef mutable struct Process 
    N = 100
    xrange::Vector{Float64} = [8500, 11500] # plot range in x (nm)
    yrange::Vector{Float64} = [5000, 15000] # plot range in y (nm) 
end

function process()
    proc = Process()

    ind = Vector{Int}(undef, 0)
    msd = Vector{Float64}(undef, proc.N-1)
    pos = Matrix{Float64}(undef, 0, 2)
    Nc = 0

    for n in 1:proc.N
        name = "./map/cornea_map_20by20_" * lpad(n, 3, '0') * ".dat"
        df = CSV.read(name, DataFrame; header=false)
        xpos = df.Column1 .* 1e9
        ypos = df.Column2 .* 1e9
        if n == 1
            for m in eachindex(xpos)
                if xpos[m] > proc.xrange[1] && xpos[m] < proc.xrange[2] &&
                    ypos[m] > proc.yrange[1] && ypos[m] < proc.yrange[2]
                    push!(ind, m)
                    # push!(pos, [xpos[m], ypos[m]])
                    pos = vcat(pos, [xpos[m] ypos[m]])
                end
            end
            Nc = length(ind)
        else
            msd[n-1] = sum((xpos[ind].-pos[:,1]).^2 + (ypos[ind].-pos[:,2]).^2)/Nc
        end
        
    end
    
    return msd
end

msd = process()
save("msd.jld", "msd", msd)
