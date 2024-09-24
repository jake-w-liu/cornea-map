using BatchAssign
using LinearAlgebra
using Infiltrator
using PlotlyJS
using FFMPEG

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
    xrange = [0, 2000] # plot range in x (nm)
    yrange = [0, 2000] # plot range in y (nm) 
end

mutable struct Cornea
    pos::Vector{Float64}
end

mutable struct Map
    grid::Matrix{Int}
end

function generate_pos(sx, sy, d)
    return [(sx - d) * rand() + par.Dc / 2, (sy - d) * rand() + d / 2]
end

function init_cornea(par::Parameters)
    CorneaList = Vector{Cornea}(undef, 0)
    pass = false
    for n = 1:par.Nc
        pos = generate_pos(par.Sx, par.Sy, par.Dc)
        if n == 1
            push!(CorneaList, Cornea(pos))
        else
            while !pass
                for c in CorneaList
                    if norm(c.pos - pos) <= par.r
                        pos .= generate_pos(par.Sx, par.Sy, par.Dc)
                        break
                    end
                    if c == CorneaList[end]
                        pass = true
                    end
                end
            end
            push!(CorneaList, Cornea(pos))
            pass = false
        end
        println(n)
    end
    return CorneaList
end

function update_cornea!(CorneaList::Vector{Cornea}, d::Real)
    for c in CorneaList
        @all c.pos[1] c.pos[2] += (rand()-0.5)*d
    end
    return nothing
end

function layout_cornea(CorneaList::Vector{Cornea}, par::Parameters)
    layout = Layout(
        template = "plotly_white",
        plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
        # xaxis=attr(scaleanchor="y"),
        xaxis = attr(range = par.xrange, title_text = "x (nm)"),
        yaxis = attr(range = par.yrange, title_text = "y (nm)"),
        shapes = [
        ],
        height = diff(par.xrange)[1]/5,
        width  = diff(par.yrange)[1]/5,
    )
    layout_include!(layout, CorneaList, par)
    return layout
end

function layout_include!(layout::PlotlyJS.Layout, CorneaList::Vector{Cornea}, par::Parameters)
    layout.shapes = []
    for n in 1:par.Nc
        if CorneaList[n].pos[1] > par.xrange[1] &&
            CorneaList[n].pos[1] < par.xrange[2] &&
            CorneaList[n].pos[2] > par.yrange[1] &&
            CorneaList[n].pos[2] < par.yrange[2]
            push!(
                layout.shapes,
                circle(
                    xref = "x",
                    yref = "y",
                    fillcolor = "rgb(251,226,81)", # KIHADA
                    line_width = 0,
                    x0 = CorneaList[n].pos[1] - par.Dc / 2,
                    y0 = CorneaList[n].pos[2] - par.Dc / 2,
                    x1 = CorneaList[n].pos[1] + par.Dc / 2,
                    y1 = CorneaList[n].pos[2] + par.Dc / 2,
                ),
            )
        end
    end
    return nothing
end

function update_plot!(plt::PlotlyJS.SyncPlot, CorneaList::Vector{Cornea}, par::Parameters)
    layout_include!(plt.plot.layout, CorneaList, par)
    react!(plt, plt.plot.data, plt.plot.layout)
end

function plot_cornea(CorneaList::Vector{Cornea}, par::Parameters)
    layout = layout_cornea(CorneaList, par)
    plt = plot(scatter(), layout)
    return plt
end

## script

# par = Parameters()
# @time CorneaList = init_cornea(par)

# layout = layout_cornea(CorneaList, par)
# plt = plot(scatter(), layout)
plt = plot_cornea(CorneaList, par)
display(plt)

if !isdir("tmp")
    mkdir("tmp")
end

for nt in 1:100
    update_cornea!(CorneaList, 20)
    update_plot!(plt, CorneaList, par)
    filename = "./tmp/snap_" * lpad(nt, 3, '0') * ".png"
    savefig(plt, filename; 
        height = round(Int, diff(par.xrange)[1]/5),
        width  = round(Int, diff(par.yrange)[1]/5),)
    # sleep(0.1)
end

imagesdirectory = "./tmp"
framerate = 30
gifname = "output.gif"
FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -i $(imagesdirectory)/snap_%03d.png -y $(gifname)`)