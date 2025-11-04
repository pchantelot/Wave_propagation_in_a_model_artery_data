using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using FileIO, UnPack, NaturalSort
cd(@__DIR__)
include("../preamble.jl")

# Material properties
ρf = 1e3
g = 9.81
sc = 40e-2 / 2090
# Height Field
H = load(joinpath(@__DIR__,"heightfield/selectedfield.jld2"))
@unpack h = H
# Line fields
p2names = filter(x -> occursin("jld2",x), readdir(joinpath(@__DIR__,"lines/")))
p2names = p2names[[2,4,1,3]]
p2color_list = [1, 5, 1, 5]
# Dispersion relation
p3names = filter(x -> occursin("jld2",x), readdir(joinpath(@__DIR__,"dispersion/")))
p3names = reverse(p3names)
p3color_list = [5, 1]

with_theme(My_theme, palette = (color = ColorSchemes.Blues_7, marker = [:circle])) do 
    fig = Figure(size = (2 * 402, 2 / 3 * 402), figure_padding = (16,16,1,1)) 
        
        # Panel 1: Schematic
        panel1 = fig[1,1] = GridLayout()
        Label(panel1[1, 1, TopLeft()], "(a)", padding = (-88, 0, -10, 0), fontsize = 16)
            ax = Axis(panel1[1,1])
            ax.limits = (-31, 66, -9, 25)
            ax.aspect = DataAspect()
            hidedecorations!(ax)
            hidespines!(ax)
                # waveguide
                # reservoir
                poly!(ax, Point2f[(50,4.5),(53,4.5),(53,13),(52,13),(52,5.5),(50,5.5)], color = :skyblue1, strokecolor = :black, strokewidth = 2)
                poly!(ax, Point2f[(50.5,13),(54.5, 13),(54.5,21),(50.5,21)], color = :skyblue1, strokecolor = :black, strokewidth = 2)
                poly!(ax, Point2f[(50,19.5),(55, 19.5),(55,21.5),(50,21.5)], color = :white)
                # right 
                x = 0:0.1:50    
                y = 10 .+ 3*sin.(x .- 5.5) .* exp.(-0.08* (x .- 5.5)) .* (x .> 5.5)
                band!(ax, x, fill(0,length(x)), y; color = :skyblue1)
                lines!(ax, Point2f[(0.3,10),(0.3,0),(49.7,0),(49.7,10)], color = :black, linewidth = 3)
                lines!(ax, x, y, color = :hotpink2, linewidth = 4.5)
                #poly!(ax, Rect(-0.3, 9.5, 50.6, 1), color = :hotpink2)
                text!(ax, 49, 11, text = "Air", align = (:right, :bottom), fontsize = 16)
                text!(ax, 49, 9, text = "Water", align = (:right, :top), fontsize = 16)
                # left
                poly!(ax, Rect(-17, 0, 10, 10), color = :skyblue1, 
                    strokecolor = :black, strokewidth = 3)
                z = -5.3:0.1:5.3
                tb = 0.1
                w = 10.6
                y =  130*(1 .- 4/(w*tb)*coth(tb*w/2) .- (2*z/w).^2 .+ 4/(w*tb)*cosh.(tb*z)/sinh(tb*w/2))
                band!(ax, z .- 12, fill(9.5,length(z)), 10 .+ y; color = :skyblue1)
                lines!(ax, z .- 12, 10 .+ y, color = :hotpink2, linewidth = 5)
                # deflection
                arrows!(ax, [-12], [9.1], [0], [3], linewidth = 2, arrowsize = 10)
                text!(ax, -12.5, 16, text = L"\delta(z)", align = (:right, :center), fontsize = 16)
                # right arrow
                arrows!(ax, [57], [12], [0], [6.5], linewidth = 2, arrowsize = 12)
                arrows!(ax, [57], [18], [0], [-7.5], linewidth = 2, arrowsize = 12)
                text!(ax, 58, 15, text = L"\Delta H", align = (:left, :center), fontsize = 16)
                # left arrow
                arrows!(ax, [-19.5], [2], [0], [7], linewidth = 2, arrowsize = 12)
                arrows!(ax, [-19.5], [8], [0], [-7], linewidth = 2, arrowsize = 12)
                text!(ax, -20.7, 5.4, text = L"$1$ cm", align = (:right, :center), fontsize = 16)
                arrows!(ax, [-5.2], [12], [0], [-2], linewidth = 2, arrowsize = 12)
                arrows!(ax, [-5.2], [8], [0], [2], linewidth = 2, arrowsize = 12)
                text!(ax, -4.5, 10, text = L"h_0", align = (:left, :center), fontsize = 16)
                # camera
                poly!(ax, Rect(23,20, 4, 4), color = :black)
                text!(ax, 28, 22, text = "Camera", align = (:left, :center), fontsize = 16)
                lines!(ax, [22, 23.5], [18, 20.5], color = :black, linewidth = 1.5)
                lines!(ax, [26.5, 28], [20.5, 18], color = :black, linewidth = 1.5)
                # gravity
                arrows!(ax, [13], [24], [0], [-4], linewidth = 2, arrowsize = 12)
                text!(ax, 14.5, 23, text = L"\mathbf{g}", align = (:left, :center), fontsize = 16)
                # Bottom arrows
                lines!(ax,[0, 50], [-1.5, -1.5], color = :black, linestyle = (:dot,:dense))
                arrows!(ax, [5], [-3.5], [44], [0], linewidth = 2, arrowsize = 12)
                arrows!(ax, [45], [-3.5], [-44], [0], linewidth = 2, arrowsize = 12)
                text!(ax, 25, -6, text = L"$58$ cm", align = (:center, :center), fontsize = 16)
                arrows!(ax, [-12], [-2.5], [4], [0], linewidth = 2, arrowsize = 12)
                arrows!(ax, [-12], [-2.5], [-4], [0], linewidth = 2, arrowsize = 12)
                text!(ax, -12, -6, text = L"$w_0 = 1$ cm", align = (:center, :center), fontsize = 16)
                # shaker
                arrows!(ax, [4.5], [10], [0], [3], linewidth = 2, arrowsize = 12)
                arrows!(ax, [4.5], [10], [0], [-3], linewidth = 2, arrowsize = 12)
                text!(ax, 4.5, 16, text = "Shaker", align = (:center, :center), fontsize = 16)
                # displacement
                arrows!(ax, [13.25], [8], [0], [3], linewidth = 2, arrowsize = 10)
                text!(ax, 14, 11.5, text = L"u_y(x,z,t)", align = (:left, :bottom), fontsize = 16)
                # Axis 
                arrows!(ax, [3.5, 3.5], [2, 2], [0, 4], [4, 0], linewidth = 1.5, arrowsize = 10)
                text!(ax, 9.7, 2, text = L"x", align = (:right, :center), fontsize = 16)
                text!(ax, 3.2, 8, text = L"y", align = (:right, :center), fontsize = 16)
                arrows!(ax, [-9.5, -9.5], [2, 2], [0, -4], [4, 0], linewidth = 1.5, arrowsize = 10)
                text!(ax, -16, 2, text = L"z", align = (:left, :center), fontsize = 16)
                text!(ax, -10, 8.5, text = L"y", align = (:left, :center), fontsize = 16)

        # Panel 2: Height Field
        panel2 = fig[2,1] = GridLayout()
        Label(panel2[1, 1, TopLeft()], "(b)", padding = (-48, 0, 7, 0), fontsize = 16)
            ax21 = Axis(panel2[1,1])
            #ax21.alignmode = Outside()
            ax21.ylabel = L"$u_y$ (arb. u.)"
            ax21.xlabel = L"x \,\, \mathrm{(cm)}"
            ax21.limits = (-1, 20, -1, 8.7)
            #ax21.yticks  = [-1,0,1]
            for i in eachindex(p2names)
                data = load(joinpath(@__DIR__,"lines/",p2names[i]))
                @unpack line = data
                    lines!(ax21, (0:length(line)-1)*sc*1e2, 2.5*(i-1) .+ line, 
                        color = reverse(ColorSchemes.Blues_7)[p2color_list[i]])
            end
            ax22 = Axis(panel2[1,2])
            ax22.limits = (-1, 20, -1, 8.7)
            hidedecorations!(ax22)
            hidespines!(ax22)
            colsize!(panel2, 2, Relative(1/12))
            colgap!(panel2, 0)
            bracket!(ax22, 1, 0, 1, 2.5,
                text = "50 Hz", textoffset = 5, width = 10, orientation = :down, 
                color = :black, textcolor = :black, fontsize = 14)
            bracket!(ax22, 1, 7.5, 1, 5,
                text = "150 Hz", textoffset = 5, width = 10, orientation = :up, 
                color = :black, textcolor = :black, fontsize = 14)

        panel3 = fig[1:2,2] = GridLayout()
        Label(panel3[1, 1, TopLeft()], "(c)", padding = (-35, 0, -10, 0), fontsize = 16)
            ax3 = Axis(panel3[1,1])
            #ax3.alignmode = Outside()
            ax3.limits = (0, 1000, 0, 350)
            ax3.ylabel = L"$f$ (Hz)"
            ax3.xlabel = L"k_x \,\, \mathrm{(1/m)}"
            for i in eachindex(p3names)[1:2]
                data = load(joinpath(@__DIR__,"dispersion/",p3names[i]))
                @unpack f, k, ΔP = data
                scatter!(ax3, k[f .> 10], f[f .> 10], 
                    markersize = 14, color = reverse(ColorSchemes.Blues_7)[p3color_list[i]],
                    label = L"\Delta P = %$(Int(round(ρf*g*ΔP*1e-2; sigdigits = 2))) \,\,\mathrm{Pa}")
            end
            axislegend(ax3, position = :lt, labelsize = 16)
    rowsize!(fig.layout, 1, Relative(3/5))
    rowgap!(fig.layout, 0)
    display(fig)
    save(joinpath(@__DIR__,"Figure2.pdf"),fig; pt_per_unit = 1)
end  