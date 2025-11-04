using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using FileIO, UnPack, NaturalSort
cd(@__DIR__)
include("../preamble.jl")

cfiles = filter(x -> occursin("disp",x),readdir(joinpath(@__DIR__,"courbure/")))
cfiles = sort(cfiles, lt = natural)
cpredic = filter(x -> occursin("_2",x),readdir(joinpath(@__DIR__,"courbure/")))
cpredic = sort(cpredic, lt = natural)

with_theme(My_theme, palette = (color = reverse(ColorSchemes.RdPu_4), marker = [:circle])) do 
    fig = Figure(size = (402, 2/3*402), figure_padding = (2,14,1,10))

        ax22 = Axis(fig[1,1])
        ax22.ylabel = L"f \,\, \mathrm{(Hz)}"
        ax22.xlabel = L"k_x \,\, \mathrm{(1/m)}"
        ax22.limits = (0, 800, 0, 350)
        for i in eachindex(cfiles) 
            data = load(joinpath(@__DIR__,"courbure/",cfiles[i]))
            @unpack f, k = data
            scatter!(ax22, k[real(k) .> 80], f[real(k) .> 80], markersize = 14, color = ColorSchemes.YlGn_4[i+1])
        end
        start = reverse([200,100,20])
        lab = reverse(["0.7", "1.0", "∞"])
        for i in eachindex(cpredic) 
            data = load(joinpath(@__DIR__,"courbure/",reverse(cpredic)[i]))
            @unpack f,k,s = data
            lines!(ax22, k, f, color = tuple.(reverse(ColorSchemes.YlGn_4)[i], s ./ maximum(s)), label = L"R/w_0 = %$(lab[i])" => (; color = reverse(ColorSchemes.YlGn_4)[i]))
        end
        axislegend(ax22, position = :rb, labelsize = 16)

        ax21 = Axis(fig[1,1])
        ax21.aspect = DataAspect()
        ax21.halign = -0.2
        ax21.valign = 1
        ax21.width = Relative(3/5)
        ax21.height = Relative(3/5)
        ax21.limits = (-17, 32, -20, 25)
        hidedecorations!(ax21)
        hidespines!(ax21)
        # Arc
        arc!(ax21,Point2f(0),15, -3*pi/2 , pi/2 - 1/2, linewidth = 2 , color = :black, linestyle = :dash)
        arc!(ax21,Point2f(0), 15, pi/2 -1/2, pi/2 + 1/2, linewidth = 10 , color = :hotpink2)
        arc!(ax21,Point2f(0), 15, pi/2 + 1/2, pi/2 + 1/2 + 3/10, linewidth = 10 , color = :grey25)
        arc!(ax21,Point2f(0), 15, pi/2 - 1/2 - 3/10, pi/2 -1/2, linewidth = 10 , color = :grey25)
        arc!(ax21,Point2f(0), 19, pi/2 -1/2 + 1/10, pi/2 + 1/2 - 1/10, color = :black, linewidth = 2)
        arrows!(ax21, 19*[cos(pi/2-1/2+1/10)], 19*[sin(pi/2-1/2+1/10)], 
            2*[cos(1/2-1/10)], 2*[-sin(1/2-1/10)], arrowsize = 12)
        arrows!(ax21, 19*[cos(pi/2+1/2-1/10)], 19*[sin(pi/2+1/2-1/10)], 
            2*[-cos(-1/2+1/10)], 2*[sin(-1/2+1/10)], arrowsize = 12)
        text!(ax21, 0, 20, text = L"w_0", align = (:center, :bottom), fontsize = 16)
        # Radius
        lines!(ax21, [0, 15], [0, 0] , color = :black)
        text!(ax21, 10, -1, text = L"R", align = (:center, :top), fontsize = 16)
        #Axis
        poly!(ax21, Circle(Point2f(16, 20), 2.5), 
            strokecolor = :black, strokewidth = 1.3, color = (:white, 0))
        poly!(ax21, Circle(Point2f(16, 20), .5), 
            strokecolor = :black, strokewidth = 1.3, color = (:black, 1))
        text!(ax21, 19, 22, text = L"\mathbf{e}_x", 
            align = (:left, :center), fontsize = 16)
        arrows!(ax21, 15 * cos(pi/16) * [1, 1], 15 * sin(pi/16) * [1, 1], 6*[cos(pi/16), -sin(pi/16)], 6*[sin(pi/16), cos(pi/16)],
            linewidth = 1.5, arrowsize = 10)
        text!(ax21, 21.5, 2.5, text = L"r", align = (:left, :center), fontsize = 16)
        text!(ax21, 14, 11, text = L"\theta", align = (:left, :center), fontsize = 16)
        
    display(fig)
    save(joinpath(@__DIR__,"Figure4.pdf"),fig; pt_per_unit = 1)
end      
