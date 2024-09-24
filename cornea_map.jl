using BatchAssign
using LinearAlgebra
using Infiltrator
using PlotlyJS

# basic parameters
Base.@kwdef mutable struct Parameters
    Dc = 29.74 # cornea diameter (nm)
    Nc = 86775 # number of cornea fibers
    Sx = 30000 # space lengh in x (nm)
    Sy = 22000 # space lengh in y (nm)
    ds = 10 # grid size (nm)
    r = 59.48 # spacing between cornea fibers (nm)
    Nx = round(Int, Sx / ds)
    Ny = round(Int, Sy / ds)
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

# initialize cornea fibers

function init_cornea(par::Parameters)
    CorneaList = Vector{Cornea}(undef, 0)
    pass = false
    for n in 1:par.Nc
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

# script
par = Parameters()
@time CorneaList = init_cornea(par)

## test1

layout = Layout(
    template="plotly_white",
    # plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
    # xaxis=attr(scaleanchor="y"),
    xaxis=attr(range=[0, par.Sx/1]),
    yaxis=attr(range=[0, par.Sy/1]),
    shapes=[
        rect(
            x0=0, y0=0, x1=par.Sx, y1=par.Sy,
            line=attr(
                color="rgb(12,72,66)", # ONANDO
                width=2,
            ),
            fillcolor="rgb(43,95,117)", # NOSHIMEHANA
            xref='x', yref='y'
        )
    ],
    height=par.Sy/20,
    width=par.Sx/20,
)
for n in 1:par.Nc
    push!(layout.shapes, circle(xref="x", yref="y",
        fillcolor="rgb(251,226,81)", # KIHADA
        # line_color="rgb(255,196,8)", # TOHOH
        x0=CorneaList[n].pos[1] -par.Dc/2, y0=CorneaList[n].pos[2] -par.Dc/2, 
        x1=CorneaList[n].pos[1] +par.Dc/2, y1=CorneaList[n].pos[2] +par.Dc/2,)
    )
end

plt = plot(layout)
display(plt)

## test2

layout = Layout(
    template="plotly_white",
    # plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
    # xaxis=attr(scaleanchor="y"),
    xaxis=attr(range=[0, par.Sx]),
    yaxis=attr(range=[0, par.Sy]),
    shapes=[
        rect(
            x0=0, y0=0, x1=par.Sx, y1=par.Sy,
            line=attr(
                color="rgb(12,72,66)", # ONANDO
                width=2,
            ),
            fillcolor="rgb(43,95,117)", # NOSHIMEHANA
            xref='x', yref='y'
        )
    ],
    height=par.Sy/par.Dc,
    width=par.Sx/par.Dc,
)
@all x y = zeros(par.Nc)
for n in 1:par.Nc
    x[n] = CorneaList[n].pos[1]
    y[n] = CorneaList[n].pos[2]
end
trace = scatter(
    x=x,
    y=y,
    mode="markers",
    marker=attr(
        color="rgb(251,226,81)", # KIHADA
        size=10,
    )
)

plt = plot(trace, layout)
