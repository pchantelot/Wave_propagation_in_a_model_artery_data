using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using FileIO, UnPack, NaturalSort
cd(@__DIR__)
include("../preamble.jl")

with_theme(My_theme) do 
    fig = Figure(size = (502, 220), figure_padding = (1,14,1,10))
        ax21 = Axis(fig[1,1])
        ax21.aspect = DataAspect()
        ax21.limits = (-12, 12, -2.5, 16)
        hidedecorations!(ax21)
        hidespines!(ax21)
        # Define deformed profile
        R = 12
        θp = pi/2 - 2/3
        θm = pi/2 + 2/3
        z = -R*cos(θp):0.01:R*cos(θp)
        tb = 0.1
        w = 2*R*cos(θp)
        y =  80*(1 .- 4/(w*tb)*coth(tb*w/2) .- (2*z/w).^2 .+ 4/(w*tb)*cosh.(tb*z)/sinh(tb*w/2)) .+ R .- z.^2/(2*R)
        yref = R .- z.^2/(2*R)
        band!(ax21, z, fill(0,length(z)), y; color = (:skyblue1, 0.5))
        band!(ax21, z, fill(0,length(z)), yref; color = (:skyblue1, 1))
        # Axis
        arrows!(ax21, [R*cos(θp)+2], [0], [-19], [0], linewidth = 3, arrowsize = 16)
        lines!(ax21, [R*cos(θp), R*cos(θp)], 0.5*[-1, 1], color = :black)
        lines!(ax21, -[R*cos(θp), R*cos(θp)], 0.5*[-1, 1], color = :black)
        arrows!(ax21, [0], [0], [0], [7], linewidth = 3, arrowsize = 16)
        text!(ax21, R*cos(θp), -1.5, text = L"-w'/2", align = (:center, :center), fontsize = 16)
        text!(ax21, -R*cos(θp), -1.5, text = L"w'/2", align = (:center, :center), fontsize = 16)
        text!(ax21, -11, -0.5, text = L"z", align = (:center, :center), fontsize = 20)
        text!(ax21, 1, 8, text = L"y", align = (:center, :center), fontsize = 20)
        # Cone def
        lines!(ax21, [0, R*cos(θp)], [0, R*sin(θp)], linewidth = 2, color = :black, linestyle = :dash)
        lines!(ax21, [0, -R*cos(θp)], [0, R*sin(θp)], linewidth = 2, color = :black, linestyle = :dash)
        arc!(ax21,Point2f(0), R, θp, θm, linewidth = 6 , color = :hotpink2)
        
        lines!(ax21, z, y, linewidth = 6 , color = (:hotpink2, 0.75))
        
        # Angles
        arc!(ax21,Point2f(0), 4, θp, pi/2, linewidth = 2 , color = :black)
        text!(ax21, 1.7, 5, text = L"\beta", align = (:center, :center), fontsize = 20)
        # Radii
        text!(ax21, 5.5, 6, text = L"\bar{R}", align = (:center, :top), fontsize = 20)
        # Pressure
        text!(ax21, 0, 10.5, text = L"P^\infty", align = (:center, :center), fontsize = 20)
        text!(ax21, -0.35, 13.1, text = L"P^\infty + P'", align = (:center, :center), fontsize = 20)

        
    display(fig)
    save(joinpath(@__DIR__,"FigureC.pdf"),fig; pt_per_unit = 1)
end      

## Saving pie syntax for the future
# Reference and deformed pies
#pie!(ax21, [0.8] , normalize=false, offset = pi/2 - 1/2 + 1/10, radius = rdr, color = (:skyblue1, 0.5), strokewidth = 0)
#pie!(ax21, [1] , normalize=false, offset = pi/2 - 1/2, radius = 15, color = :skyblue1)