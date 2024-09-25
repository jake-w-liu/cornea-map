using BatchAssign
using LinearAlgebra
using Infiltrator
using PlotlyJS
using FFMPEG
using Base.Threads

# basic parameters
Base.@kwdef mutable struct Parameters
    Dc = 29.74 # cornea diameter (nm)
    Nc = 86775 # number of cornea fibers
    Sx = 22000 # space lengh in x (nm)
    Sy = 30000 # space lengh in y (nm)
    ds = 10 # grid size (nm)
    r = 59.48 # spacing between cornea fibers (nm)
    Nx = round(Int, Sx / ds)
    Ny = round(Int, Sy / ds)
    xrange = [0, 500] # plot range in x (nm)
    yrange = [0, 500] # plot range in y (nm) 
end

mutable struct CorneaList
    pos::Matrix{Float64}
    nb::Vector{Vector{Int}}
end

function generate_pos(sx, sy, d)
    return [(sx - d) * rand() + par.Dc / 2, (sy - d) * rand() + d / 2]
end

function find_neighbor(pos::Matrix{<:Real})
    N = size(pos, 1)
    neighbor = Vector{Vector{Int}}(undef, N)
    dif = zeros(N)
    si = zeros(Int, N)
    
    for nc in 1:N
        @threads for n in 1:N
            dif[n] = norm(pos[nc] .- pos[n])
        end
        si .= sortperm(dif)
        neighbor[nc] = [si[2], si[3], si[4]]
        println(nc)
    end
    return neighbor
end

function init_cornea(par::Parameters)
    cor_pos = zeros(par.Nc, 2)
    pass = false
    for n = 1:par.Nc
        pos = generate_pos(par.Sx, par.Sy, par.Dc)
        if n == 1
            cor_pos[n, :] .= pos
        else
            while !pass
                for nc in 1:n
                    if norm(cor_pos[nc] .- pos) <= par.r
                        pos .= generate_pos(par.Sx, par.Sy, par.Dc)
                        break
                    end
                    if nc == n
                        pass = true
                    end
                end
            end
            cor_pos[n, :] .= pos
            pass = false
        end
        println(n)
    end
    
    cor_ind = find_neighbor(cor_pos)
    cl = CorneaList(cor_pos, cor_ind)
    return cl
end

function update_cornea!(cl::CorneaList, par::Parameters, d::Real)
    pass = false
    tmp = zeros(2)
    for n in 1:size(cl.pos, 1)
        pass = false
        while !pass
            # @all cl.pos[n, 1] cl.pos[n, 2] += (rand()-0.5)*d
            tmp[1] = cl.pos[n, 1] + (rand()-0.5)*d
            tmp[2] = cl.pos[n, 2] + (rand()-0.5)*d
            if norm(tmp.- cl.pos[cl.nb[n][1], :]) > par.r &&  
                norm(tmp .- cl.pos[cl.nb[n][2], :]) > par.r &&  
                norm(tmp .- cl.pos[cl.nb[n][3], :]) > par.r 
                cl.pos[n, :] .= tmp
                pass = true
            end
        end
    end
    return nothing
end

function layout_cornea(cl::CorneaList, par::Parameters)
    layout = Layout(
        template = "plotly_white",
        plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
        # xaxis=attr(scaleanchor="y"),
        xaxis = attr(range = par.xrange, title_text = "x (nm)"),
        yaxis = attr(range = par.yrange, title_text = "y (nm)"),
        shapes = [
        ],
        height = 500,
        width  = round(Int, 500/diff(par.xrange)[1]*diff(par.yrange)[1]),
    )
    layout_include!(layout, cl, par)
    return layout
end

function layout_include!(layout::PlotlyJS.Layout, cl::CorneaList, par::Parameters)
    layout.shapes = []
    for n in 1:par.Nc
        if cl.pos[n, 1] > par.xrange[1] &&
            cl.pos[n, 1] < par.xrange[2] &&
            cl.pos[n, 2] > par.yrange[1] &&
            cl.pos[n, 2] < par.yrange[2]
            push!(
                layout.shapes,
                circle(
                    xref = "x",
                    yref = "y",
                    fillcolor = "rgb(251,226,81)", # KIHADA
                    line_width = 0,
                    x0 = cl.pos[n, 1] - par.Dc / 2,
                    y0 = cl.pos[n, 2] - par.Dc / 2,
                    x1 = cl.pos[n, 1] + par.Dc / 2,
                    y1 = cl.pos[n, 2] + par.Dc / 2,
                ),
            )
        end
    end
    return nothing
end

function update_plot!(plt::PlotlyJS.SyncPlot, cl::CorneaList, par::Parameters)
    layout_include!(plt.plot.layout, cl, par)
    react!(plt, plt.plot.data, plt.plot.layout)
end

function plot_cornea(cl::CorneaList, par::Parameters)
    layout = layout_cornea(cl, par)
    plt = plot(scatter(), layout)
    return plt
end

## script

# par = Parameters()
# @time cl = init_cornea(par)

plt = plot_cornea(cl, par)
display(plt)

if !isdir("tmp")
    mkdir("tmp")
end

for nt in 1:100
    update_cornea!(cl, par, 20)
    update_plot!(plt, cl, par)
    filename = "./tmp/snap_" * lpad(nt, 3, '0') * ".png"
    savefig(plt, filename; 
        height = 500,
        width  = round(Int, 500/diff(par.xrange)[1]*diff(par.yrange)[1]),)
    # sleep(0.1)
    println(nt)
end

framerate = 10
gifname = "output.gif"
FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -i ./tmp/snap_%03d.png -y $(gifname)`)

# if isdir("tmp")
#     rm("tmp"; recursive=true)
# end