using CairoMakie, ColorSchemes, LaTeXStrings, MathTeXEngine
using FileIO, UnPack, NaturalSort
using Statistics
using Integrals
using Roots
using Dierckx
cd(@__DIR__)
include("../preamble.jl")

# Parameters (Ecoflex OO-30)
E = 3 * 23e3
h0 = 0.77e-3
w0 = 1.e-2
g = 9.81
ρf = 1e3
# water height range
x = 0.01:0.1:50

# DL solution for initial guess
P = ρf*g*x*1e-2
T = (4/3*E*h0*w0^2*P.^2/24).^(1/3)
δ = P*w0^2 ./ (8*T)
λ3 = 1 .+ P.^2*w0^2 ./(24*T.^2)
κ =  P ./ T #This is the maximum curvature (i.e. second derivative)

# Reference computation: linear elasticity (no P = λP)
δb = Vector{Float64}()
λ3b = Vector{Float64}()
κbm = Vector{Float64}()
Tb = Vector{Float64}()
Ab = Vector{Float64}()

for i  in eachindex(x)

    Ptemp = P[i]
    λ30 = λ3[i] # initial guess based on DL solution

    function myintegral(lbd) # computes ℓ/w

        a =  3 * w0 / (2*h0) * sqrt(4/3*(lbd-1))
        toint(z,p) = 1/2 * sqrt(1 + (3*Ptemp*w0/(8*E*h0)/(lbd-1)*(-z+sinh(a*z)/sinh(a)))^2)

        result = solve(IntegralProblem(toint, (-1.0, 1.0)), QuadGKJL()).u 
    
        return(result)
    
    end

    f = u -> u - myintegral(u)
    λtemp = fzero(f, λ30)
    Ttemp = 4/3* E * h0 * (λtemp-1)
    atemp = 3 * w0 / (2*h0) * sqrt(4/3*(λtemp-1))

    # Compute average curvature 
    bd = 1/atemp * acosh(sinh(atemp)/atemp)
    toint(y,p) = 1/(2*bd) * (Ptemp/Ttemp .* (1 .- atemp*cosh.(atemp*y)/sinh(atemp))) ./
    (1 .+ (Ptemp/Ttemp*w0/2 .*(.-y .+ sinh.(atemp*y)./sinh(atemp))).^2).^(3/2)
    result = solve(IntegralProblem(toint, (-bd, bd)), QuadGKJL()).u

    # Compute area 
    tointA(y,p) = w0/2*(Ptemp*w0^2)/(8*Ttemp) * (1 -  2/atemp * coth(atemp) - 
        y^2 + 2/atemp * cosh(atemp*y) / sinh(atemp))
    resultA = solve(IntegralProblem(tointA, (-1, 1)), QuadGKJL()).u

    push!(Tb, Ttemp)
    push!(λ3b, λtemp)
    push!(δb, Ptemp*w0^2/(8*Ttemp)*(1+2/atemp*(1-cosh(atemp))/sinh(atemp)))
    push!(κbm,result)
    push!(Ab,resultA)

end

with_theme(My_theme) do 
    f = Figure(size = (402, 176), figure_padding = (3,14,1,3))
    panel1 = f[1,1] = GridLayout()
        ax21 = Axis(panel1[1,1])
        ax21.xlabel = L"\Delta P \,\, \mathrm{(Pa)}"
        ax21.ylabel = L"\lambda"
        ax21.limits = (0, 3200, 1, 1.17)
            lines!(ax21, P, λ3b, color = :black)
            #lines!(ax21, P, λ3b)
        ax22 = Axis(panel1[1,2])
        ax22.xlabel = L"\Delta P \,\, \mathrm{(Pa)}"
        ax22.ylabel = L"\bar{R}/w_0"
        ax22.limits = (0, 3200, 0.5, 1.2)
            lines!(ax22, P, 1. ./ (κbm*w0), color = :black)

    display(f)
    save(joinpath(@__DIR__,"FigureB.pdf"), f)
end

# Saving the associated Moens Korteweg velocity
# 1. linear elasticity
splA = Spline1D(x, 1 .+ Ab/w0^2)
dAdP = derivative(splA, x; nu = 1)
vMKlin = sqrt.((w0^2 .+ Ab) ./ (ρf * w0^2/(ρf*g*1e-2) .* dAdP))
save(joinpath(@__DIR__,"vMKlin.jld2"),"ΔP",P,"v",vMKlin)
# 2. equivalent strip model
vb = Vector{Float64}()
for i in eachindex(x)

    Rs =  1 ./ κbm[i]
    λs = λ3b[i]
    ws = 2 * Rs .* sin.(λs * w0 ./ (2*Rs))
    As = Ab[i]
    hs = h0 / λs
    Ts = 4/3*E*hs*(λs-1)
    # Compute reference length in approx geometry
    toint0(z,p) = ws / 2 * sqrt(1 + (z*ws/(2*Rs))^2)
    l0 = solve(IntegralProblem(toint0, (-1.0, 1.0)), QuadGKJL()).u 

    Ac = Vector{Float64}()
    for j in eachindex(x[1:20])

        Ptmp = P[j]

        function myintegral(lbd) # computes ℓ/w

            a =  3 * ws / (2 * hs) * sqrt(4/3 * (λs-1))
            toint(z,p) = ws / 2 * sqrt((1 + (ws/2 * (Ptmp / Ts - (4/3 * E * hs * (lbd-1)) / (Rs*Ts))*(-z+sinh(a*z)/sinh(a)) -z*ws/(2*Rs))^2))
            result = solve(IntegralProblem(toint, (-1.0, 1.0)), QuadGKJL()).u 
        
            return(result)
        
        end

        f = u -> u - myintegral(u)/l0
        λtemp = fzero(f, 1 + 3 / (4*E*hs)*Ptmp*Rs - 0.0001)
        Ttemp = 4/3 * E * hs * (λtemp-1)
        atemp = 3 * ws / (2*hs) * sqrt(4/3 * (λs-1))

        # Compute area 
        tointA(y,p) = ws/2*(ws^2)/8 * (Ptmp/Ts - Ttemp/(Rs*Ts))* (1 -  2/atemp * coth(atemp) - 
            y^2 + 2/atemp * cosh(atemp*y) / sinh(atemp))
        resultA = solve(IntegralProblem(tointA, (-1, 1)), QuadGKJL()).u
        push!(Ac, resultA)

    end

    splAR = Spline1D(x[1:20], 1 .+ Ac/w0^2)
    spldAdPR = derivative(splAR, x[1]; nu = 1)*w0^2/(ρf*g*1e-2)
    push!(vb, sqrt.((w0^2 .+ As) ./ (ρf * spldAdPR)))

end

save(joinpath(@__DIR__,"vMKdef.jld2"),"ΔP",P,"v",vb)
