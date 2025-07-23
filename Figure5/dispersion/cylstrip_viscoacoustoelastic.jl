# Guided wave dispersion in a stretched viscoelastic strip in cylindrical coordinates
# Load packages
# Easy identity matrix and kronecker tensor notation
using LinearAlgebra, Kronecker
# SCM 
push!(LOAD_PATH, string(pwd(),raw"/Documents/GitHub/JuliaSCM/"))
using SCM
# Plotting 
using CairoMakie

## Problem definition: We look for Lamb waves propagating in direction 1 in a curved strip. 
# 
# Strip with height h, width w and curvature radius R in vacuum
#
#                           top, zero normal traction
#  e2 (= r)             ---------------------------------- 
#  |_ e3 (= θ)     left - solid: μs, λs, ρs              - right, zero displacement
#                       ---------------------------------- 
#                           bottom, zero normal traction
#
# Waves propagate in direction 1 (= x),  U(x,t) = u(r,θ)e^{i(kx-ωt)}
# normalization using ρs, μs and h
# Choose degrees of freedom
udof = 1:3 # displacement dof, coupled Lamb and SH waves

# Solid properties (Ecoflex OO-30)
ρs = 1070
μs = 23e3
λs = 1000^2*ρs - 2*μs # measured speed of sound in Ecoflex 1000 m/s
h = 0.77e-3
w = 1e-2
R = 0.78*w # radius of curvature
# Define parameters to compute the 4th order elastic tensor
λ1 = 1.
λ3 = 1.06
τ = 330e-6
n = 0.32
α = 0.29
β = 0.29
# update geometry (λ1λ2λ3 = 1)
hdef = h / (λ1*λ3)
wdef = w*λ3

# Discretization: We differentiate in direction 2 and 3using SCM.
N = 5 # number of nodes in direction 2, must be even for in-plane modes and odd for out-of-plane modes
P = 10 # number of nodes in direction 3
x2, Dt1 = chebdif(N,2) # Nodes in [-1, 1]
x2 = (R .+ (hdef .* x2) ./ 2) ./ hdef
D2 = 2 .* (hdef/hdef) * Dt1[:,:,1]
D22 = 4 .* (hdef/hdef)^2 * Dt1[:,:,2]
_, Dt2 = chebdif(P,2)
D3 = 2 * (R/wdef) .* Dt2[:,:,1] 
D33 = 4 *(R/wdef)^2 .* Dt2[:,:,2]
A = [0 0 0; 0 0 -1; 0 1 0] # Derivation of the basis vectors
D2d = I(P) ⊗ D2; 
D22d = I(P) ⊗ D22; 
x2dinv = I(P) ⊗ Diagonal(1 ./ x2); 
x2dinv2 = I(P) ⊗ Diagonal(1 ./ x2.^2); 
D3d = D3 ⊗ I(N); 
D33d = D33 ⊗ I(N); 

## Solve polynomial eigenvalue problem
# Impose ω and get eigenvalues and eigenvectors
ω = 2 * pi .* range(50, 400, length = 60) 
k = Array{ComplexF64,2}(undef,length(ω),2*length(udof)*N*P)
u = Array{ComplexF64,3}(undef,length(ω),length(udof)*N*P,2*length(udof)*N*P)

for i in axes(ω,1)

    C = C_viscoacoustoelastic(λ1, λ3, sqrt(μs/ρs), sqrt((λs+2*μs)/ρs), ρs, ω[i] / (2*pi), τ, n, α, β)
    C = C ./ μs

    # Material matrices
    C11 = C[1,udof,udof,1]
    C12 = C[1,udof,udof,2]
    C21 = C[2,udof,udof,1]
    C13 = C[1,udof,udof,3]
    C31 = C[3,udof,udof,1]
    C22 = C[2,udof,udof,2]
    C23 = C[2,udof,udof,3]
    C32 = C[3,udof,udof,2]
    C33 = C[3,udof,udof,3]
    M = collect(I(length(udof)))

    # Discretized matrices to build the polynomial eigenvalue problem
    # (ik)^2 Lkk + (ik) Lk + L0 + ω^2 Md = 0
    Lkk = C11 ⊗ I(N*P); Lkk = collect(Lkk)
    Lk = (C12 .+ C21) ⊗ D2d .+ C21 ⊗ x2dinv .+ 
        C13*I(3) ⊗ (x2dinv*D3d) .+ C13*A ⊗ x2dinv .+ 
        I(3)*C31 ⊗ (x2dinv*D3d) .+ A*C31 ⊗ x2dinv 
    L0 = C22 ⊗ D22d .+ C22 ⊗ (x2dinv*D2d) .+ 
        C23*I(3) ⊗ (D2d*x2dinv*D3d) .+ C23*A ⊗ (D2d*x2dinv) .+ 
        I(3)*C32 ⊗ (D3d*x2dinv*D2d) .+ A*C32 ⊗ (x2dinv*D2d) .+ 
        C23*I(3) ⊗ (x2dinv2*D3d) .+ C23*A ⊗ x2dinv2 .+
        I(3)*C33*I(3) ⊗ (x2dinv2*D33d) .+ I(3)*C33*A ⊗ (x2dinv2*D3d) .+ 
        A*C33*I(3) ⊗ (x2dinv2*D3d) .+ A*C33*A ⊗ x2dinv2
    Md = M ⊗ I(N*P); Md = collect(Md)

    # boundary conditions
    # N×P vector to N×P matrix
    bc = reshape(1:N*P,N,P)
    top = bc[1,:]
    bottom = bc[N,:]
    left = bc[:,1]
    right = bc[:,P]
    # boundary traction operators top
    Bkt = C21 ⊗ I(N*P)[top,:]
    B0t = C22 ⊗ D2d[top,:] .+ C23*I(3) ⊗ (x2dinv*D3d)[top,:] .+ C23*A ⊗ x2dinv[top,:]
    # boundary traction operators bottom
    Bkb = C21 ⊗ I(N*P)[bottom,:]
    B0b = C22 ⊗ D2d[bottom,:] .+ C23*I(3) ⊗ (x2dinv*D3d)[bottom,:] .+ C23*A ⊗ x2dinv[bottom,:]
    # zero displacement left
    B0l = I(length(udof)) ⊗ I(N*P)[left,:]
    # zero displacement right
    B0r = I(length(udof)) ⊗ I(N*P)[right,:]

    # Impose boundary conditions
    topall = (0:length(udof)-1)*N*P .+ top'; topall = topall'[:]
    bottomall = (0:length(udof)-1)*N*P .+ bottom'; bottomall = bottomall'[:]
    leftall = (0:length(udof)-1)*N*P .+ left'; leftall = leftall'[:]
    rightall = (0:length(udof)-1)*N*P .+ right'; rightall = rightall'[:]

    # top
    Lkk[topall,:] .= 0
    Lk[topall,:] .= Bkt
    L0[topall,:] .= B0t
    Md[topall,:] .= 0
    # bottom
    Lkk[bottomall,:] .= 0
    Lk[bottomall,:] .= Bkb
    L0[bottomall,:] .= B0b
    Md[bottomall,:] .= 0
    # left
    Lkk[leftall,:] .= 0
    Lk[leftall,:] .= 0
    L0[leftall,:] .= B0l
    Md[leftall,:] .= 0
    # right
    Lkk[rightall,:] .= 0
    Lk[rightall,:] .= 0
    L0[rightall,:] .= B0r
    Md[rightall,:] .= 0

    local temp = ω[i] * hdef / sqrt(μs/ρs)
    local problem = PEP([temp^2 .* Md .+ L0, Lk, Lkk])
    local ktmp, utmp = polyeig(problem)
    k[i,:] = -1im * ktmp' / hdef
    u[i,:,:] = utmp

end

# Store result of a 2D calculation before (and for) postprocessing
result = result2D(udof,N,P,ω,k,u)
idx = modedirection(result)

# remove negative real(k) and purely imaginary solutions
f = repeat(ω,1,Int(2*(length(udof)*N*P))) / (2*pi)
k[(0.49*pi .< angle.(k)) .| (angle.(k) .< -0.49*pi)] .= NaN 
k[imag(k) .> 0 ] .= NaN
k[real(k) .> 10e2] .= NaN
k[real(k) .< 1e2] .= NaN
k[idx .!== 2] .= NaN
f[isnan.(k)] .= NaN
idx[isnan.(k)] .= 4

## Plots
kplot = k[:]
fplot = f[:]
filter!(!isnan,kplot)
filter!(!isnan,fplot)
shading  = 1 .- abs.(imag.(kplot)) / (0.05*maximum(abs.(imag.(kplot))))

# extra filtering only for display
shading[shading .< 0] .= NaN
fplot[isnan.(shading)] .= NaN
kplot[isnan.(shading)] .= NaN
filter!(!isnan,fplot)
filter!(!isnan,kplot)
filter!(!isnan,shading)
# Sort
p = sortperm(fplot)
ff = fplot[p]
kk = kplot[p]
ss = shading[p]
# Keep only first modes
A = unique(ff)
knew = []
snew = []
for i in eachindex(A)
    idxtemp = findall(x-> x== A[i],ff)
    maxtemp, maxidx = findmax(real(kk[idxtemp]))
    push!(knew,maxtemp)
    push!(snew,ss[idxtemp[maxidx]])
end

with_theme(theme_latexfonts()) do 

    fig = Figure(size = (1.2 * 360, 360))

        ax = Axis(fig[1,1])
        ax.limits = (0, 800, 0, 350)
        ax.xlabel = L"k \,\, \mathrm{(1/m)}"
        ax.ylabel = L"f \,\, \mathrm{(Hz)}"
        ax.xlabelsize = 20 
        ax.ylabelsize = 20
        ax.title = "Out of plane modes"
        scatter!(ax,real(kplot),fplot, 
        color = tuple.(:red, shading))
        #scatter!(ax,real(knew),A)
        colsize!(fig.layout, 1, Aspect(1, 1.2))

    resize_to_layout!(fig)
    display(fig)
end
R = "∞"
save(joinpath(@__DIR__,"0030_8cm.jld2"),"f",A,"k",knew,"s",snew)