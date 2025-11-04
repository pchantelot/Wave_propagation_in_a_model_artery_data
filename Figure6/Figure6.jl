using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using FileIO, UnPack, NaturalSort
cd(@__DIR__)
include("../preamble.jl")

# Experimental and predicted dispersion relations
wnames = filter(x -> occursin("Data",x), readdir(joinpath(@__DIR__,"dispersion/")))
wnames = sort(wnames, lt = natural)
pnames = filter(x -> occursin("simple.jld2",x), readdir(joinpath(@__DIR__,"dispersion/")))
sort!(pnames, lt = natural)
# Measured Vb
dP = [0, 2, 5, 8, 15, 30] *1000*9.81*1e-2
Vb = [2.8, 3.5, 3.8, 4.1, 4.8, 5.1]

with_theme(My_theme, palette = (color = reverse(ColorSchemes.Blues_7), marker = [:circle])) do 
    fig = Figure(size = (4/2 *402, 2/3 * 402), figure_padding = (4,16,3,1))
        
        panel2 = fig[1,1] = GridLayout()
        Label(panel2[1, 1, TopLeft()], "(a)", padding = (-40, 0, -5, 0), fontsize = 16)
        Label(panel2[1, 2, TopLeft()], "(b)", padding = (-40, 0, -5, 0), fontsize = 16)
            ax21 = Axis(panel2[1,1])
            ax21.ylabel = L"f \,\, \mathrm{(Hz)}"
            ax21.xlabel = L"k_x \,\, \mathrm{(1/m)}"
            ax21.limits = (0, 1000, 0, 350)
            for name in wnames[1:1]
                data = load(joinpath(@__DIR__,"dispersion/",name))
                @unpack f, k = data
                scatter!(ax21, k[f .> 10][1:2:end], f[f .> 10][1:2:end], markersize = 14)
            end
            for name in pnames[1:1]
                data = load(joinpath(@__DIR__,"dispersion/",name))
                @unpack f,k = data
                lines!(ax21, real(k)[k .> 100], f[k .> 100])
            end
            # Inset
            inset = Axis(panel2[1,1])
            inset.halign = 0.24
            inset.valign = 0.93
            inset.width = Relative(45/100)
            inset.height = Relative(45/100)
            inset.limits = (0, 150, 0, 50)
            inset.xticklabelsize = 12
            inset.xticks = [0, 100]
            inset.yticklabelsize = 12
            inset.yticks = [0, 25, 50]
            #hidedecorations!(inset, ticks = false)
            x = 0:1000
            data = load(joinpath(@__DIR__,"dispersion/",wnames[1]))
            @unpack f, k = data
                scatter!(inset, k[f .> 10], f[f .> 10], markersize = 14)
                lines!(inset, x, 2.8*x/(2*pi), linestyle = :dash)
                poly!(inset, Point2f[(45, 25), (70, 25+2.8*25/(2*pi)), (45, 25+2.8*25/(2*pi))],
                    color = (:grey, 0), strokewidth = 2, strokecolor = :black)
                text!(inset, -1, 43, text = L"V_{b}", align = (:left, :center), 
                    fontsize = 20, color = :black)
            translate!(inset.scene, 0, 0, 10)
            # this needs separate translation as well, since it's drawn in the parent scene
            translate!(inset.elements[:background], 0, 0, 9)

            ax22 = Axis(panel2[1,2])
            ax22.ylabel = L"V_{b} \,\, \mathrm{(m/s)}"
            ax22.xlabel = L"Δ P \,\, \mathrm{(Pa)}"
            ax22.limits = (-150, 3200, 0, 8)
            data = load(joinpath(@__DIR__,"vMKlin.jld2"))
            @unpack v, ΔP = data
                lines!(ax22, ΔP, v, color = :black, linestyle = :dot, label = "static model")
            data = load(joinpath(@__DIR__,"vMKdef.jld2"))
            @unpack v, ΔP = data
                    lines!(ax22, ΔP, v, color = :black, label = "incremental\n model")
                scatter!(ax22, dP, Vb,
                    color = 1:6, colormap = reverse(ColorSchemes.Blues_7))
            axislegend(ax22, position = :rb)
            rowgap!(panel2, 0)

            panel1 = fig[1,2] = GridLayout()
            Label(panel1[1, 1, TopLeft()], "(c)", padding = (-40, 0, -5, 0), fontsize = 16)
                ax1 = Axis(panel1[1,1])
                ax1.ylabel = L"f \,\, \mathrm{(Hz)}"
                ax1.xlabel = L"k_x \,\, \mathrm{(1/m)}"
                ax1.limits = (0, 1000, 0, 350)
                ax1.xticks = ([0, 100, 500, 1000], ["0", L"1/w_0", "500", "1000"])
                pressures = Vector{LaTeXString}()
                poly!(ax1, Point2f[(0, 0), (100, 0), (100, 350), (0, 350)],
                    color = (:grey, 0.2), strokewidth = 0, strokecolor = :black)
                #text!(ax1, 500, 325, text = L"k_x w_0 >1", align = (:center, :center), 
                #    fontsize = 20, color = :black)
                for name in wnames
                    data = load(joinpath(@__DIR__,"dispersion/",name))
                    @unpack f, k, ΔP = data
                    scatter!(ax1, k[f .> 10][1:2:end], f[f .> 10][1:2:end], markersize = 14)
                    push!(pressures, L"%$(Int(round(ΔP*1000*9.81*1e-2, sigdigits = 2)))")
                end
                for name in pnames
                    data = load(joinpath(@__DIR__,"dispersion/",name))
                    @unpack f,k = data
                    lines!(ax1, real(k)[k .> 100]*0.93, f[k .> 100])
                end
                # group legend for the pressure values, labels are in pressures
                group_p = [MarkerElement(marker = :circle, color = reverse(ColorSchemes.Blues_7)[i], strokecolor = :black, strokewidth = 1, markersize = 16) for i in eachindex(pressures)]
                # labels for predictions
                lab = [L"1.0, \,  \infty", L"1.02, \, 1.25", 
                    L"1.04, \, 0.9", L" 1.06, \,  0.78", L" 1.09, \,  0.66", L" 1.15, \,  0.58"]
                group_l = [LineElement(linestyle = :solid, color = reverse(ColorSchemes.Blues_7)[i], linewidth = 3) for i in eachindex(lab)]
                axislegend(ax1, group_l, lab, L"\lambda_\theta, \,\, \bar{R}/w_0:"; position = (1.01, -0.1), 
                    labelsize = 12, orientation = :horizontal, nbanks = 3, patchsize = (15, 20), 
                    titlesize = 18, titlehalign= :right, titlegap = 1, rowgap = 0)
    colsize!(fig.layout, 2, Relative(1/2))
    display(fig)
    save(joinpath(@__DIR__,"Figure6.pdf"),fig; pt_per_unit = 1)
end      

