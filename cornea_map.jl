using BatchAssign
using LinearAlgebra
using Infiltrator
using PlotlyJS
using FFMPEG
using Base.Threads
using BenchmarkTools
using Printf
using JLD

# basic parameters
Base.@kwdef mutable struct Parameters
    Dc::Float64 = 30 # cornea diameter (nm)
    ri::Float64 = 1.504 # refractive index
    # rho::Float64 = 0.00013147727272727272 # density
    rho::Float64 = 0.000000013147727272727272 # density ## Test
    Sx::Float64 = 20000 # space lengh in x (nm)
    Sy::Float64 = 20000 # space lengh in y (nm)
    St::Float64 = 1000 # source thickness (nm)
    Sw::Float64 = 4000 # source width (nm)
    Sp::Vector{Float64} = [500, 10000] # source center position (nm)
    Cb::Float64 = 100 # min distance from cylinder to boundary (nm)
    Nc::Int = round(Int, rho * Sx * Sy) # number of cornea fibers
    ds::Float64 = 10 # grid size (nm)
    r::Float64 = 60 # spacing between cornea fibers (nm)
    Nx::Int = round(Int, Sx / ds)
    Ny::Int = round(Int, Sy / ds)
    xrange::Vector{Float64} = [8500, 11500] # plot range in x (nm)
    yrange::Vector{Float64} = [5000, 15000] # plot range in y (nm) 
    drift::Float64 = 20 # drift motion 
    mn::Int = 10 # map number
end

mutable struct CorneaList
    pos::Matrix{Float64}
    nb::Vector{Vector{Int}}
    pos_ini::Matrix{Float64}
end

function find_neighbor(pos::Matrix{<:Real}, par::Parameters)
    N = size(pos, 1)
    neighbor = Vector{Vector{Int}}(undef, N)
    rd = zeros(2)

    @inbounds @views for n = 1:N
        neighbor[n] = []
        for nc = 1:N
            if nc != n
                rd .= pos[nc] .- pos[n]
                if norm(rd) < (par.r + 2 * par.drift)
                    pushfirst!(neighbor[n], nc)
                end
            end
        end
        println(n)
    end
    return neighbor
end

function gen_position!(par::Parameters, pos::Vector, check::Bool = false)
    check = false
    while !check
        pos[1] = (par.Sx - 2 * par.Cb) * rand() + par.Cb
        pos[2] = (par.Sy - 2 * par.Cb) * rand() + par.Cb
        if pos[1] > par.Sp[1] - par.St / 2 - par.Cb &&
           pos[1] < par.Sp[1] + par.St / 2 + par.Cb &&
           pos[2] > par.Sp[2] - par.Sw / 2 - par.Cb &&
           pos[2] < par.Sp[2] + par.Sw / 2 + par.Cb
            check = false
        else
            check = true
        end
    end

    return nothing
end

function init_cornea(par::Parameters)
    cor_pos = zeros(par.Nc, 2)
    @all pass check = false
    @all rd pos = zeros(2)
    @inbounds @views for n = 1:par.Nc
        gen_position!(par, pos, check)
        if n == 1
            cor_pos[n, :] .= pos
        else
            while !pass
                for nc = 1:n
                    rd .= cor_pos[nc, :] .- pos
                    if norm(rd) <= par.r
                        gen_position!(par, pos, check)
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

    cor_ind = find_neighbor(cor_pos, par)
    cl = CorneaList(cor_pos, cor_ind, copy(cor_pos))
    return cl
end

function update_cornea!(cl::CorneaList, par::Parameters)
    pass = false
    @all rd tmp = zeros(2)
    @inbounds @views for n in 1:par.Nc
        pass = false
        while !pass
            tmp[1] = cl.pos_ini[n, 1] + (rand() - 0.5) * par.drift * 2
            tmp[2] = cl.pos_ini[n, 2] + (rand() - 0.5) * par.drift * 2
            if tmp[1] < par.Cb &&
               tmp[1] > par.Sx - par.Cb &&
               tmp[2] < par.Cb &&
               tmp[2] > par.Sy - par.Cb
                continue
            end

            if !isempty(cl.nb[n]) 
                for nn = 1:length(cl.nb[n])
                    rd .= tmp .- cl.pos[cl.nb[n][nn], :]
                    if norm(rd) < par.r
                        break
                    end
                    if nn == length(cl.nb[n])
                        cl.pos[n, :] .= tmp
                        pass = true
                    end
                end
            else
                pass = true
            end
        end
    end
    return nothing
end

function layout_cornea(cl::CorneaList, par::Parameters)
    layout = Layout(
        template = "plotly_white",
        # plot_bgcolor = "rgb(43,95,117)", # NOSHIMEHANA
        # xaxis=attr(scaleanchor="y"),
        xaxis = attr(range = par.xrange, title_text = "x (nm)", scaleanchor="y", zeroline=false),
        yaxis = attr(range = par.yrange, title_text = "y (nm)", zeroline=false),
        shapes = [
            rect(
            x0=par.xrange[1], y0=par.yrange[1], x1=par.xrange[2], y1=par.yrange[2],
            line=attr(
                color="RoyalBlue",
                width=1,
            ),
            fillcolor="rgb(43,95,117)",
            xref='x', yref='y'
        )
        ],
        height = 1000,
        width = round(Int, 1000 / diff(par.yrange)[1] * diff(par.xrange)[1] * 1.2),
    )
    layout_include!(layout, cl, par)
    return layout
end

function layout_include!(layout::PlotlyJS.Layout, cl::CorneaList, par::Parameters)
    layout.shapes = [rect(
        x0=par.xrange[1], y0=par.yrange[1], x1=par.xrange[2], y1=par.yrange[2],
        line=attr(
            color="RoyalBlue",
            width=1,
        ),
        fillcolor="rgb(43,95,117)",
        xref='x', yref='y'
    )]
    @inbounds @views for n = 1:par.Nc
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
    return nothing
end

function plot_cornea(cl::CorneaList, par::Parameters)
    layout = layout_cornea(cl, par)
    plt = plot(scatter(), layout)
    return plt
end

function output_data(
    cl::CorneaList,
    par::Parameters,
    filename::String,
    data::Vector{Vector{Float64}},
)
    # data = Vector{Vector{Float64}}(undef, par.Nc)
    formatted_row = Vector{String}(undef, 4)
    @views for n in 1:par.Nc
        data[n] = [cl.pos[n, 1] * 1e-9, cl.pos[n, 2] * 1e-9, par.Dc / 2 * 1e-9, par.ri]
    end
    open(filename, "w") do file
        @views for row in data
            formatted_row .= [@sprintf("%e", x) for x in row]
            println(file, join(formatted_row, ","))
        end
    end
end

function output_grid(
    cl::CorneaList,
    par::Parameters,
    filename::String,
    bg_ind::Int,
    mtr_ind::Int,
    grid::Matrix,
)
    # grid = ones(Int, par.Nx, par.Ny) .* bg_ind
    fill!(grid, bg_ind)
    x = range(0, (par.Nx - 1) * par.ds, par.Nx) .+ par.ds / 2
    y = range(0, (par.Ny - 1) * par.ds, par.Ny) .+ par.ds / 2
    @all indx indy = 0
    nr = ceil(Int, par.Dc / par.ds / 2 + 1)
    @views for n = 1:par.Nc
        indx = round(Int, (cl.pos[n, 1] - x[1]) / par.ds) + 1
        indy = round(Int, (cl.pos[n, 2] - y[1]) / par.ds) + 1
        for nx = -nr:nr
            for ny = -nr:nr
                if indx + nx > 0 &&
                   indx + nx <= par.Nx &&
                   indy + ny > 0 &&
                   indy + ny <= par.Ny
                    if (x[indx+nx] - cl.pos[n, 1])^2 + (y[indy+ny] - cl.pos[n, 2])^2 <
                       (par.Dc / 2)^2
                        grid[indx+nx, indy+ny] = mtr_ind
                    end
                end
            end
        end
    end
    save(filename, "grid", grid)
end

## script

function main()
    par = Parameters()
    @time cl = init_cornea(par)

    plt = plot_cornea(cl, par)
    display(plt)

    if !isdir("tmp")
        mkdir("tmp")
    end

    if !isdir("map")
        mkdir("map")
    end

    if !isdir("map_jld")
        mkdir("map_jld")
    end

    file_name = "cornea_map_$(round(Int, par.Sx/1000))by$(round(Int, par.Sy/1000))_"
    data = Vector{Vector{Float64}}(undef, par.Nc)
    grid = ones(Int, par.Nx, par.Ny) .* 3
    for nt = 1:par.mn
        @time update_cornea!(cl, par)
        update_plot!(plt, cl, par)

        savefig(
            plt,
            "./tmp/snap_" * lpad(nt, 3, '0') * ".png";
            height = 1000,
            width = round(Int, 1000 / diff(par.yrange)[1] * diff(par.xrange)[1] * 1.2),
        )

        @time output_data(cl, par, "./map/" * file_name * lpad(nt, 3, '0') * ".dat", data)
        @time output_grid(
            cl,
            par,
            "./map_jld/" * file_name * lpad(nt, 3, '0') * ".jld",
            3,
            8,
            grid,
        )

        println(nt)
    end

    framerate = 10
    gifname = "output.gif"
    FFMPEG.ffmpeg_exe(
        `-framerate $(framerate) -f image2 -i ./tmp/snap_%03d.png -y $(gifname)`,
    )
    return nothing
end

main()



# if isdir("tmp")
#     rm("tmp"; recursive=true)
# end
