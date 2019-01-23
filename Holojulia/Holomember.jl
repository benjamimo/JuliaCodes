include("/mnt/bendata/Documents/Docs/PhysicsLocal/Beams/Beams.jl")  # include Small-core and LG beams
include("/mnt/bendata/Documents/Docs/PhysicsLocal/Beams/MatrixTrans.jl")

using Base.Threads

"""
     memberPC(Useed, Xs, Ys, w0, TopologicalCharge, RadialOrder, circle, N)

Simulates a member of the ensamble to generate a partially coherent Useed beam """
function memberPC(Useed, ren, col, circle::Float64, N::Int64)

# Preallocates memory for matrices U and Uaux
# ren = size(Ys,1)
# col = size(Xs,2)
Uaux = zeros(ComplexF64, ren, col)
U = zeros(ComplexF64, ren, col)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = Integer.(div.(r.*cos.(theta),16E-6))
yv = Integer.(div.(r.*sin.(theta),16E-6))

# Coherent superposition
for member=1:N
    Uaux = exp(im * phi[member]) .* MatrixTranslate(Useed, xv[member], yv[member])
    Uaux = Uaux./maximum(abs.(Uaux))
    U = Uaux + U
end

return U
end

"""
     memberPCmt(Useed, Xs, Ys, w0, TopologicalCharge, RadialOrder, circle, N)

Simulates a member of the ensamble to generate a partially coherent Useed beam (multithreading)"""
function memberPCmt(Useed, ren, col, circle::Float64, N::Int64)

# Preallocates memory for matrices U and Uaux
# ren = size(Ys,1)
# col = size(Xs,2)
Uaux = zeros(ComplexF64, ren, col)
U = zeros(ComplexF64, ren, col)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = Integer.(div.(r.*cos.(theta),16E-6))
yv = Integer.(div.(r.*sin.(theta),16E-6))

# Coherent superposition
@Threads.threads for member=1:N
    Uaux = exp(im * phi[member]) .* MatrixTranslate(Useed, xv[member], yv[member])
    Uaux = Uaux./maximum(abs.(Uaux))
    U = Uaux + U
end

return U
end
