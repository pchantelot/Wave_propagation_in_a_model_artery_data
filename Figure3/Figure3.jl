using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine, MakieTeX
using FileIO, UnPack, NaturalSort
cd(@__DIR__)
include("../preamble.jl")

tfiles = filter(x -> occursin("disp",x),readdir(joinpath(@__DIR__,"tension/")))
tpredic = filter(x -> occursin("_2",x),readdir(joinpath(@__DIR__,"tension/")))

with_theme(My_theme, palette = (color = reverse(ColorSchemes.Reds_4), marker = [:circle])) do 
    fig = Figure(size = (502, 502), figure_padding = (3,14,0,-10))
        panel1 = fig[1,1] = GridLayout()
        Label(panel1[1, 1, TopLeft()], "(a)", padding = (-90, 0, -40, 0), fontsize = 16)
            ax11 = Axis(panel1[1,1])
            ax11.limits = (-47, 72, -5, 35)
            ax11.aspect = DataAspect()
            hidedecorations!(ax11)
            hidespines!(ax11)
            # First view
            poly!(ax11, Rect(0, 0, 50, 10*1.3), color = :hotpink2)
            poly!(ax11, Rect(0, -4, 50, 4), color = :grey25)
            poly!(ax11, Rect(0, 10*1.3, 50, 4), color = :grey25)
            poly!(ax11, Rect(0, 1.5, 50, 10), color = (:hotpink2, 0),  strokewidth = 3, 
                strokecolor = :black, linestyle = (:dot, :dense))
            arrows!(ax11, [53], [6], [0], [5], linewidth = 2, arrowsize = 12)
            arrows!(ax11, [53], [6], [0], [-5], linewidth = 2, arrowsize = 12)
            text!(ax11, 54, 6.5, text = L"w = \lambda_z w_0", align = (:left, :center), fontsize = 16)
            #text!(ax11, 54, 13, text = L"\mathbf{F} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \lambda_3^{-1} & 0 \\ 0 & 0 & \lambda_3 \end{pmatrix}", 
            #    align = (:left, :center), fontsize = 16)
            teximg!(ax11, L"\mathbf{F} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1/\lambda_z & 0 \\ 0 & 0 & \lambda_z \end{pmatrix}";
                position = (42, 26), align = (:left, :center), scale = 1)
            # axis right
            arrows!(ax11, [0, 0], [20, 20], [0, 4], [4, 0], linewidth = 1.5, arrowsize = 10)
            text!(ax11, 4.5, 19, text = L"x", align = (:left, :center), fontsize = 16)
            text!(ax11, -1, 25.5, text = L"z", align = (:right, :center), fontsize = 16)
            # second view
            poly!(ax11, Rect(-26, 5, 10*1.3, 3), color = :hotpink2)
            poly!(ax11, Rect(-24.5, 4.55, 10, 3*1.3), color = (:hotpink2,0), 
                strokewidth = 3, strokecolor = :black, linestyle = (:dot, :dense))
            poly!(ax11, Rect(-30, 5, 4, 3), color = :grey25)
            poly!(ax11, Rect(-13, 5, 4, 3), color = :grey25)
            arrows!(ax11, [-31], [3.5], [0], [0.2], linewidth = 2, arrowsize = 10)
            arrows!(ax11, [-31], [9.5], [0], [-0.2], linewidth = 2, arrowsize = 10)
            text!(ax11, -35, 6.5, text = L"h = \frac{h_0}{λ_z} ", align = (:right, :center), fontsize = 16)
            # axis top
            arrows!(ax11, [-7, -7], [3, 3], [0, -4], [4, 0], linewidth = 1.5, arrowsize = 10)
            text!(ax11, -6, 9, text = L"y", align = (:left, :center), fontsize = 16)
            text!(ax11, -12, 2, text = L"z", align = (:right, :center), fontsize = 16)
            # Laser sheet
            θ = -30 *pi /180
            c = [-37, 15.5]
            rot = [cos(θ) -sin(θ); sin(θ)  cos(θ)]
            lines!(ax11, [-19.5 ,-31],[8, 15], color = :lime)
            poly!(Point2f[c] .+ Point2f[rot*([-41, 16] .- c), rot*([-31, 16] .- c), rot*([-31, 20] .- c ), rot*([-41, 20] .- c)], color = :black)
            x = 0:0.1:50 
            lines!(ax11, x, 6.5 .+ 4*cos.(0.5*(x .- x[1])) .* exp.(- 0.07 .* (x .- x[1])), color = :lime)
            text!(ax11, -41, 24, text = L"\mathrm{Laser \, sheet}", align = (:left, :center), fontsize = 16)
            # camera
            poly!(ax11, Rect(23, 24, 4, 4), color = :black)
            text!(ax11, 25, 31, text = L"\mathrm{Camera}", align = (:center, :center), fontsize = 16)
            lines!(ax11, [22, 23.5], [22, 24.5], color = :black, linewidth = 1.5)
            lines!(ax11, [26.5, 28], [24.5, 22], color = :black, linewidth = 1.5)

        panel2 = fig[2,1] = GridLayout()
        Label(panel2[1, 1, TopLeft()], "(b)", padding = (-35, 0, -5, 0), fontsize = 16)
        ax12 = Axis(panel2[1,1])
        #ax12.alignmode = Outside()
        ax12.ylabel = L"f \,\, \mathrm{(Hz)}"
        ax12.xlabel = L"k_x \,\, \mathrm{(1/m)}"
        ax12.limits = (0, 800, 0, 350)
        for name in tfiles  
            data = load(joinpath(@__DIR__,"tension/",name))
            @unpack f, k = data
            scatter!(ax12, k, f, markersize = 14)
        end
        for i in eachindex(tpredic)
            data = load(joinpath(@__DIR__,"tension/",tpredic[i]))
            @unpack f,k,s = data
            lines!(ax12, k, f, color = tuple.(reverse(ColorSchemes.Reds_4)[i], s ./ maximum(s)),
                label = L"\lambda_z = %$(parse(Float32,tpredic[i][6:end-7])/100)" => (; color = reverse(ColorSchemes.Reds_4)[i]))
        end
        axislegend(ax12, position = :rb, labelsize = 16)
    
        rowsize!(fig.layout, 2, Relative(602/1000))
    rowgap!(fig.layout, 0)
    display(fig)
    save(joinpath(@__DIR__,"Figure3.pdf"),fig; pt_per_unit = 1)
end      
