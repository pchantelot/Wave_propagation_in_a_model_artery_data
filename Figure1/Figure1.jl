using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using FileIO, UnPack, NaturalSort
cd(@__DIR__)
include("../preamble.jl")

# Material properties
μs = 23e3
ρs = 1070
ρf = 1000
h = 1e-3
R = 4e-3 

panel1files = filter(x -> occursin("tube",x), readdir(joinpath(@__DIR__,"dispersion/")))
panel2files = filter(x -> occursin("tube",x), readdir(joinpath(@__DIR__,"zoom_dispersion/")))
# Sort for colors
panel1files = panel1files[[4,2,5,1,3]]
panel2files = panel2files[[4,2,5,1,3]]

with_theme(My_theme, palette = (color = reverse(ColorSchemes.Set1_5)[[1,3,2,4,5]], marker = [:circle])) do 
    x = collect(1e-3:1100)
    fig = Figure(size = (2*502, 2/3*502), figure_padding = (2,14,2,2)) 
        
        # Panel 1: dispersion relation
        panel1 = fig[1,1] = GridLayout()
        Label(panel1[1, 1, TopLeft()], "(a)", padding = (0, -20, -5, 0), fontsize = 16)
            ax = Axis(panel1[1,1])
            ax.alignmode = Outside()
            ax.limits = (-10, 1000, 0, 300)
            ax.ylabel = L"f \,\, \mathrm{(Hz)}"
            ax.xlabel = L"k_x \,\, \mathrm{(1/m)}"
            ax.xticks = ([0, 250, 500, 1000], ["0", L"1/R", "500", "1000"])
            poly!(ax,Point2f[(0, 0), (150, 0), (150, 50), (0, 50)], 
                    color = (:grey, 0.2), strokecolor = :black, strokewidth = 1)
                for name in panel1files
                    data = load(joinpath(@__DIR__,"dispersion/",name))
                    @unpack f, k = data
                    p = sortperm(f)
                    scatter!(ax, real(k[p][1:end]), f[p][1:end], markersize = 14)
                end
                lines!(ax, x , x .* sqrt(3*μs/ρs) / (2*pi), color = :black, linewidth = 4)
                lines!(ax, x , x .* sqrt(μs/ρs) / (2*pi), color = :black, linestyle = (:dash, :dense), linewidth = 4)
                lines!(ax, x , x .* 1500 / (2*pi), color = :black, linestyle = :dashdot, linewidth = 4)
                
            # Sketch in inset
            inset = Axis(panel1[1,1])
            inset.halign = 1.25
            inset.valign = 0.12
            inset.aspect = DataAspect()
            inset.width = Relative(3/5)
            inset.height = Relative(3/5)
            hidedecorations!(inset)
            hidespines!(inset)
                # Make the outer gradient
                grad =  ColorScheme(range(colorant"skyblue1", colorant"white", length=25))
                grad = resample_cmap(grad, 25; alpha = Tuple([1. ,0.]))
                X = -80:80
                Y = -80:80
                r = [sqrt(x^2+y^2) for x in X, y in Y]
                r[r .> 80] .= 80
                r[r .< 45] .= 45
                col = 45 .- r
                # Define h arrow coordinate
                ra = [26, 54]
                θa = pi/2 * 1.3
                xa = ra .* cos.(θa)
                ya = ra .* sin.(θa)
                ura = [1, -1]
                uxa = ura .* cos.(θa)
                uya = ura .* sin.(θa)
                # Outer gradient
                #surface!(inset,X,Y,Z, color = col, colormap = reverse(grad), shading = NoShading)
                heatmap!(inset, X, Y, col, colormap = reverse(grad))
                # Circles
                poly!(inset, Circle(Point2f(0, 0), 45), color = :hotpink2)
                poly!(inset, Circle(Point2f(0, 0), 35), color = :skyblue1)
                # radius with text
                lines!(inset,[0, 40],[0,0], color = :black, linewidth = 2)
                text!(inset, 15, -12 , text = L"R", align = (:center, :center), fontsize = 20)
                # Arrows with text
                arrows!(inset, xa, ya, uxa, uya, arrowsize = 16)
                text!(inset, -10, 60 , text = L"h", align = (:center, :center), fontsize = 20)
                # x axis symbol
                poly!(inset, Circle(Point2f(-65, -40), 7), 
                    strokecolor = :black, strokewidth = 1.3, color = (:white, 0))
                lines!(inset, -65 .+ [-7, 7] .* cos(pi/4), -40 .+ [-7, 7 ] .* sin(pi/4), 
                    color = :black, linewidth = 1.3)
                lines!(inset, -65 .+ [-7, 7] .* cos(pi/4), -40 .+ [7, -7 ] .* sin(pi/4), 
                    color = :black, linewidth = 1.3)
                text!(inset, -46, -48, text = L"\mathbf{e}_x", 
                    align = (:center, :center), fontsize = 16)
                
        # Panel 2: Low frequency zoom + 
        panel2 = fig[1,2] = GridLayout()
        Label(panel2[1, 1, TopLeft()], "(b)", padding = (-40, 0, -5, 0), fontsize = 16)
            ax2 = Axis(panel2[1,1])
            ax2.limits = (-3, 150, 0, 50)
            ax2.ylabel = L"f \,\, \mathrm{(Hz)}"
            ax2.xlabel = L"k_x \,\, \mathrm{(1/m)}"
                for name in panel2files
                    data = load(joinpath(@__DIR__,"zoom_dispersion/",name))
                    @unpack f, k = data
                    p = sortperm(f)
                    scatter!(ax2, real(k[p]), f[p], markersize = 14)
                end
                lines!(ax2, x , x .* 1500 / (2*pi), color = :black, linestyle = :dashdot, linewidth = 4,
                    label = L" V_\phi = V_f")
                lines!(ax2, x, x .* sqrt(3*μs/ρs) / (2*pi), color = :black, linewidth = 4,
                    label = L" V_\phi = \sqrt{E_s/\rho_s}")
                lines!(ax2, x, x  .* sqrt(μs/ρs) / (2*pi), color = :black, linewidth = 4, 
                    linestyle = (:dash, :dense),
                    label = L"V_\phi = \sqrt{\mu_s/\rho_s}")
                lines!(ax2, x, x .* sqrt(3*μs/ρf * h / R / 2) / (2*pi),
                    label = L"V_\phi = V_{MK}, \; \mathrm{Eq. \; (3)}", color = ColorSchemes.Set1_4[2])
                lines!(ax2, x, (x).^2 .* sqrt.(3*μs*R*h ./ ρf) / (2*pi),
                    color = tuple.(ColorSchemes.Set1_4[1],0.5), linestyle = :dash)
                lines!(ax2, x, (x).^2 .* sqrt.(3*μs*R*h ./ (3*ρf)) / (2*pi),
                    label = "Eq. (4)", color = ColorSchemes.Set1_4[1])
            leg = Legend(panel2[2,1], ax2, orientation = :horizontal, nbanks = 5, labelsize = 16, tellheight = true, patchsize =(35, 20))
            leg.alignmode = Mixed(bottom = 0)
        rowgap!(panel2, 10)

        # panel3: mode shapes
        panel3 = fig[1,3] = GridLayout()
        Label(panel3[1, 1, TopLeft()], "(c)", padding = (0, 0, -5, 0), fontsize = 16)
            ax31 = Axis(panel3[2,2], alignmode = Outside(3,3,0,0))
            hidedecorations!(ax31)
            hidespines!(ax31)
            ax31.aspect = DataAspect()
                img = load(joinpath(@__DIR__,"Modedisplay/Flexuralmode.png"))
                image!(ax31,rotr90(img))
                text!(ax31, 25, 525, text = L"\mathrm{Flexion}", 
                    align = (:left, :center), fontsize = 20, color = ColorSchemes.Set1_4[1])
                Box(panel3[2,2], cornerradius = 10, color = (:white, 0), 
                    strokecolor = ColorSchemes.Set1_4[1], strokewidth = 2.6)

            ax32 = Axis(panel3[2,1], alignmode = Outside(3,3,0,0))
            hidedecorations!(ax32)
            hidespines!(ax32)
            ax32.aspect = DataAspect()
                img = load(joinpath(@__DIR__,"Modedisplay/MKmode.png"))
                image!(ax32,rotr90(img))
                text!(ax32, 25, 525, text = L"\mathrm{Breathing}", 
                    align = (:left, :center), fontsize = 20, color = ColorSchemes.Set1_4[2])
                Box(panel3[2,1], cornerradius = 10, color = (:white, 0), 
                    strokecolor = ColorSchemes.Set1_4[2], strokewidth = 2.6)

            ax33 = Axis(panel3[1,2], alignmode = Outside(3,3,0,0))
            hidedecorations!(ax33)
            hidespines!(ax33)
            ax33.aspect = DataAspect()
                img = load(joinpath(@__DIR__,"Modedisplay/Torsionmode.png"))
                image!(ax33,rotr90(img))
                text!(ax33, 25, 525, text = L"\mathrm{Torsion}", 
                    align = (:left, :center), fontsize = 20, color = ColorSchemes.Set1_4[4])
                Box(panel3[1,2], cornerradius = 10, color = (:white, 0), 
                    strokecolor = ColorSchemes.Set1_4[4], strokewidth = 2.6)

            ax34 = Axis(panel3[1,1], alignmode = Outside(3,3,0,0))
            hidedecorations!(ax34)
            hidespines!(ax34)
            ax34.aspect = DataAspect()
                img = load(joinpath(@__DIR__,"Modedisplay/Compressionmode.png"))
                image!(ax34,rotr90(img))
                text!(ax34, 25, 525, text = L"\mathrm{Compression}", 
                    align = (:left, :center), fontsize = 20, color = ColorSchemes.Set1_4[3])
                Box(panel3[1,1], cornerradius = 10, color = (:white, 0), 
                    strokecolor = ColorSchemes.Set1_4[3], strokewidth = 2.6)
        rowgap!(panel3, 12)
        colgap!(panel3, 12)

    colsize!(fig.layout, 1, Relative(2/5))
    colsize!(fig.layout, 2, Relative(1/6))
    colgap!(fig.layout, 10)
    display(fig)
    save(joinpath(@__DIR__,"Figure1.pdf"),fig; pt_per_unit = 1)
end