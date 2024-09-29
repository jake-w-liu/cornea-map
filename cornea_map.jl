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
    rho::Float64 = 0.00013147727272727272 # density
    Sx::Float64 = 22000 # space lengh in x (nm)
    Sy::Float64 = 30000 # space lengh in y (nm)
    Nc::Int = round(Int, rho * Sx * Sy) # number of cornea fibers
    ds::Float64 = 10 # grid size (nm)
    r::Float64 = 60 # spacing between cornea fibers (nm)
    Nx::Int = round(Int, Sx / ds)
    Ny::Int = round(Int, Sy / ds)
    xrange::Vector{Float64} = [0, 1000] # plot range in x (nm)
    yrange::Vector{Float64} = [0, 1000] # plot range in y (nm) 
    drift::Float64 = 20 # drift motion 
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
    # @threads for n in 1:N
    #     neighbor[n] = []
    #     @inbounds @views for nc in 1:N
    #         if nc != n
    #             # rd .= pos[nc] .- pos[n]
    #             if norm(pos[nc] .- pos[n]) < 2*(par.r + par.drift)
    #                 pushfirst!(neighbor[n], nc)
    #             end
    #         end
    #     end
    #     println(n)
    # end
    return neighbor
end

function init_cornea(par::Parameters)
    cor_pos = zeros(par.Nc, 2)
    pass = false
    @all rd pos = zeros(2)
    @inbounds @views for n = 1:par.Nc
        pos[1] = (par.Sx - par.Dc) * rand() + par.Dc / 2
        pos[2] = (par.Sy - par.Dc) * rand() + par.Dc / 2
        if n == 1
            cor_pos[n, :] .= pos
        else
            while !pass
                for nc = 1:n
                    rd .= cor_pos[nc, :] .- pos
                    if norm(rd) <= par.r
                        pos[1] = (par.Sx - par.Dc) * rand() + par.Dc / 2
                        pos[2] = (par.Sy - par.Dc) * rand() + par.Dc / 2
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
    N = size(cl.pos, 1)
    @inbounds @views for n = 1:N
        pass = false
        while !pass
            tmp[1] = cl.pos_ini[n, 1] + (rand() - 0.5) * par.drift * 2
            tmp[2] = cl.pos_ini[n, 2] + (rand() - 0.5) * par.drift * 2
            if tmp[1] < 0 && tmp[1] > par.Sx && tmp[2] < 0 && tmp[2] > par.Sy
                continue
            end

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
        shapes = [],
        height = 500,
        width = round(Int, 500 / diff(par.xrange)[1] * diff(par.yrange)[1]),
    )
    layout_include!(layout, cl, par)
    return layout
end

function layout_include!(layout::PlotlyJS.Layout, cl::CorneaList, par::Parameters)
    layout.shapes = []
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
end

function plot_cornea(cl::CorneaList, par::Parameters)
    layout = layout_cornea(cl, par)
    plt = plot(scatter(), layout)
    return plt
end

function output_data(cl::CorneaList, par::Parameters, filename::String)
    data = Vector{Vector{Float64}}(undef, par.Nc)
    formatted_row = Vector{String}(undef, 4)
    for n = 1:par.Nc
        data[n] = [cl.pos[n, 1] * 1e-9, cl.pos[n, 2] * 1e-9, par.Dc / 2 * 1e-9, par.ri]
    end
    open(filename, "w") do file
        for row in data
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
)
    grid = ones(Int, par.Nx, par.Ny) .* bg_ind
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
        # _, indx = findmin(abs.(x .-cl.pos[n, 1]))
        # _, indy = findmin(abs.(y .-cl.pos[n, 2]))
        # grid[indx, indy] = mtr_ind
        # for nx in -nr:nr
        #     for ny in -nr:nr
        #         if (nx * par.ds)^2 + (ny * par.ds)^2 < (par.Dc/2)^2
        #             if indx+nx > 0 && indx+nx <= par.Nx && 
        #                 indy+ny > 0 && indy+ny <= par.Ny
        #                 grid[indx+nx, indy+ny] = mtr_ind
        #             end
        #         end
        #     end
        # end
    end
    save(filename, "grid", grid)
end

## script

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

for nt = 1:100
    @time update_cornea!(cl, par)
    update_plot!(plt, cl, par)

    savefig(
        plt,
        "./tmp/snap_" * lpad(nt, 3, '0') * ".png";
        height = 500,
        width = round(Int, 500 / diff(par.xrange)[1] * diff(par.yrange)[1]),
    )

    @time output_data(cl, par, "./map/map_" * lpad(nt, 3, '0') * ".dat")
    @time output_grid(cl, par, "./map_jld/map_" * lpad(nt, 3, '0') * ".jld", 3, 8)

    println(nt)
end

framerate = 10
gifname = "output.gif"
FFMPEG.ffmpeg_exe(`-framerate $(framerate) -f image2 -i ./tmp/snap_%03d.png -y $(gifname)`)

# if isdir("tmp")
#     rm("tmp"; recursive=true)
# end