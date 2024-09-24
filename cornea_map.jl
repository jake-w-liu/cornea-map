using BatchAssign
using LinearAlgebra
using Infiltrator
using PlotlyJS

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
for nt in 1:100
    update_cornea!(CorneaList, 20)
    update_plot!(plt, CorneaList, par)
    sleep(0.1)
end

# ## test1

# layout = Layout(
#     template = "plotly_white",
#     # plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
#     # xaxis=attr(scaleanchor="y"),
#     xaxis = attr(range = [0, par.Sx / 1]),
#     yaxis = attr(range = [0, par.Sy / 1]),
#     shapes = [
#         rect(
#             x0 = 0,
#             y0 = 0,
#             x1 = par.Sx,
#             y1 = par.Sy,
#             line = attr(
#                 color = "rgb(12,72,66)", # ONANDO
#                 width = 2,
#             ),
#             fillcolor = "rgb(43,95,117)", # NOSHIMEHANA
#             xref = 'x',
#             yref = 'y',
#         ),
#     ],
#     height = par.Sy / 20,
#     width = par.Sx / 20,
# )
# for n = 1:par.Nc
#     push!(
#         layout.shapes,
#         circle(
#             xref = "x",
#             yref = "y",
#             fillcolor = "rgb(251,226,81)", # KIHADA
#             # line_color="rgb(255,196,8)", # TOHOH
#             x0 = CorneaList[n].pos[1] - par.Dc / 2,
#             y0 = CorneaList[n].pos[2] - par.Dc / 2,
#             x1 = CorneaList[n].pos[1] + par.Dc / 2,
#             y1 = CorneaList[n].pos[2] + par.Dc / 2,
#         ),
#     )
# end

# plt = plot(layout)
# display(plt)

# ## test2

# layout = Layout(
#     template = "plotly_white",
#     # plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
#     # xaxis=attr(scaleanchor="y"),
#     xaxis = attr(range = [0, par.Sx]),
#     yaxis = attr(range = [0, par.Sy]),
#     # shapes=[
#     #     rect(
#     #         x0=0, y0=0, x1=par.Sx, y1=par.Sy,
#     #         line=attr(
#     #             color="rgb(12,72,66)", # ONANDO
#     #             width=2,
#     #         ),
#     #         fillcolor="rgb(43,95,117)", # NOSHIMEHANA
#     #         xref='x', yref='y'
#     #     )
#     # ],
#     height = par.Sy / par.Dc,
#     width = par.Sx / par.Dc,
# )
# @all x y = zeros(par.Nc)
# for n = 1:par.Nc
#     x[n] = CorneaList[n].pos[1]
#     y[n] = CorneaList[n].pos[2]
# end
# trace = scatter(
#     x = x,
#     y = y,
#     mode = "markers",
#     marker = attr(
#         color = "rgb(251,226,81)", # KIHADA
#         size = 1,
#     ),
# )

# plt = plot(trace, layout)

# ## test3

# layout = Layout(
#     template = "plotly_white",
#     plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
#     # xaxis=attr(scaleanchor="y"),
#     xaxis = attr(range = [0, par.Sx / 1]),
#     yaxis = attr(range = [0, par.Sy / 1]),
#     shapes = [
#         circle(
#             xref = "x",
#             yref = "y",
#             fillcolor = "rgb(251,226,81)", # KIHADA
#             # line_color="rgb(255,196,8)", # TOHOH
#             x0 = [],
#             y0 = [],
#             x1 = [],
#             y1 = [],
#         ),
#     ],
#     height = par.Sy / 20,
#     width = par.Sx / 20,
# )
# for n = 1:par.Nc
#     push!(layout.shapes[1].x0, CorneaList[n].pos[1] - par.Dc / 2)
#     push!(layout.shapes[1].x1, CorneaList[n].pos[1] + par.Dc / 2)
#     push!(layout.shapes[1].y0, CorneaList[n].pos[2] - par.Dc / 2)
#     push!(layout.shapes[1].y1, CorneaList[n].pos[2] + par.Dc / 2)
# end

# plt = plot(layout)
