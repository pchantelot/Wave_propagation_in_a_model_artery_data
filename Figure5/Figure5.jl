using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using FileIO, UnPack, NaturalSort
cd(@__DIR__)
include("../preamble.jl")

afiles = filter(x -> occursin("disp",x),readdir(joinpath(@__DIR__,"dispersion/")))
apredic = filter(x -> occursin("0030",x),readdir(joinpath(@__DIR__,"dispersion/")))
profiles = filter(x -> occursin("profiles",x),readdir(joinpath(@__DIR__,"profiles/")))
w = 1e-2
lab = [L"1.0, \,  \infty", L"1.02, \, 1.25", 
    L"1.04, \, 0.9", L" 1.06, \,  0.78"]

with_theme(My_theme, palette = (color = reverse(ColorSchemes.OKeeffe2[1:2:end]), marker = [:circle])) do 
    fig = Figure(size = (502, 502), figure_padding = (3,14,1,1))
        panel1 = fig[2,1] = GridLayout()
        Label(panel1[1, 1, TopLeft()], "(b)", padding = (-40, 0, -5, 0), fontsize = 16)
            ax11 = Axis(panel1[1,1])
            ax11.ylabel = L"\delta/w_0"
            ax11.xlabel = L"2z/w_0"
            ax11.yticks = [0, 0.1, 0.2]
            ax11.limits = (-1.1, 1.1, -0.01, 0.2)
            ax11.xticks = [-1, 0, 1]
            # hidexdecorations!(ax11, ticks = false)
            ax12 = Axis(panel1[1,2])
            ax12.ylabel = L"|\kappa| w_0"
            ax12.xlabel = L"2z/w_0"
            ax12.limits = (-1.1, 1.1, -0.1, 2)
            ax12.xticks = [-1, 0, 1]
            # ΔP = 0
            lines!(ax11,[-1,1],[0,0], color = reverse(ColorSchemes.OKeeffe2[1:2:end])[1])
            # ΔP > 0
            lines!(ax12,[-1,1],[0,0], color = reverse(ColorSchemes.OKeeffe2[1:2:end])[1])
            for i in eachindex(profiles)
                data = load(joinpath(@__DIR__,"profiles/",profiles[i]))
                @unpack δ, κ, κm, z, bd = data
                lines!(ax11, z, δ/w, color = reverse(ColorSchemes.OKeeffe2[1:2:end])[i+1])
                #vspan!(ax12, -bd, bd, color = :lightgrey)
                lines!(ax12, z, κ * w, color = reverse(ColorSchemes.OKeeffe2[1:2:end])[i+1])
                lines!(ax12, bd*[-1, 1], w*κm*[1, 1], 
                    linestyle = (:dot, :dense), color = reverse(ColorSchemes.OKeeffe2[1:2:end])[i+1])
            end
            
            
        panel2 = fig[1,1] = GridLayout()
        Label(panel2[1, 1, TopLeft()], "(a)", padding = (-35, 0, -5, 0), fontsize = 16)
            ax2 = Axis(panel2[1,1])
            ax2.ylabel = L"f \,\, \mathrm{(Hz)}"
            ax2.xlabel = L"k_x \,\, \mathrm{(1/m)}"
            ax2.limits = (0, 800, 0, 350)
            for name in afiles  
                data = load(joinpath(@__DIR__,"dispersion/",name))
                @unpack f, k = data
                scatter!(ax2, k[k .> 100], f[k .> 100], markersize = 14)
            end
            for i in eachindex(apredic)
                data = load(joinpath(@__DIR__,"dispersion/",apredic[i]))
                @unpack f,k,s = data
                lines!(ax2, k, f, color = tuple.(reverse(ColorSchemes.OKeeffe2[1:2:end])[i], s ./ maximum(s)), label = lab[i] => (; color= reverse(ColorSchemes.OKeeffe2[1:2:end])[i]))
            end
            axislegend(ax2, L"λ_\theta, \,\, \bar{R}/w_0: ", position = :rb, orientation = :horizontal, nbanks = 2,
                padding = (0,0,0,0), rowgap = 0, titlehalign= :right, titlegap = 1, labelsize = 16, titlesize = 20)
            inset = Axis(panel2[1,1])
            inset.aspect = DataAspect()
            inset.width = Relative(3/5)
            inset.height = Relative(3/5)
            inset.halign = -0.25
            inset.valign = 1.15
            inset.limits = (-1, 30, -8, 25)
            hidedecorations!(inset)
            hidespines!(inset)
            # Section
            poly!(inset, Rect(0, 0, 10, 10), color = (:skyblue1,0), 
            strokecolor = :black, strokewidth = 3)
            z = -5.3:0.1:5.3
            tb = 0.1
            wd = 10.6
            y =  130*(1 .- 4/(wd*tb)*coth(tb*wd/2) .- (2*z/wd).^2 .+ 4/(wd*tb)*cosh.(tb*z)/sinh(tb*wd/2))
            band!(inset, z .+ 5, fill(9.65,length(z)), 10 .+ y; color = :white)
            lines!(inset, z .+ 5, 10 .+ y, color = :hotpink2, linewidth = 5)
            arrows!(inset, [5], [9.6], [0], [2.5], linewidth = 2, arrowsize = 10)
            text!(inset, 4.5, 15, text = L"\delta(z)", align = (:right, :center), fontsize = 16)
            # reservoir
            poly!(inset, Point2f[(15,16),(15,18),(11,18),(11,5.5),(10,5.5),(10,4.5),(12,4.5),(12,17),(14,17),(14,16)],
                color = :white, strokecolor = :black, strokewidth = 2)
            poly!(inset, Point2f[(13,11),(16,11),(16,16),(13,16)], color = :skyblue1, strokecolor = :black, strokewidth = 2)
            poly!(inset, Point2f[(13.2,12.),(15.8,12.),(15.8,15.8),(13.2,15.8)], color = :white)
            poly!(inset, Point2f[(14,4.5),(20,4.5),(20,13),(19,13),(19,5.5),(15,5.5), (15,11),(14,11),(14,4.5)], color = :skyblue1, strokecolor = :black, strokewidth = 2)
            poly!(inset, Point2f[(17.5,13),(21.5, 13),(21.5,21),(17.5,21)], color = :skyblue1, strokecolor = :black, strokewidth = 2)
            poly!(inset, Point2f[(17,19.5),(22, 19.5),(22,21.5),(17,21.5)], color = :white)
            arrows!(inset, [23], [15], [0], [3.5], linewidth = 2, arrowsize = 12)
            arrows!(inset, [23], [16.5], [0], [-3.5], linewidth = 2, arrowsize = 12)
            text!(inset, 24, 15.75, text = L"\Delta H", align = (:left, :center), fontsize = 16)
            
    
    rowsize!(fig.layout, 1, Relative(7/10))
    rowgap!(fig.layout, 0)
    display(fig)
    save(joinpath(@__DIR__,"Figure5.pdf"),fig; pt_per_unit = 1)
end      