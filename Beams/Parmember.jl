include("/mnt/bendata/Documents/Docs/PhysicsLocal/Beams/Beams.jl")  # include Small-core and LG beams
include("/mnt/bendata/Documents/Docs/PhysicsLocal/Beams/MatrixTrans.jl")

using Base.Threads

"""
     memberPC(Useed, Xs, Ys, w0, TopologicalCharge, RadialOrder, circle, N)

Simulates a member of the ensamble to generate a partially coherent Useed beam """
function memberPC(Useed, ren, col, XYsize, circle::Float64, N::Int64)

# Preallocates memory for matrices U and Uaux
Uaux = zeros(ComplexF64, ren, col)
U = zeros(ComplexF64, ren, col)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = Integer.(div.(r.*cos.(theta),XYsize))
yv = Integer.(div.(r.*sin.(theta),XYsize))

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
function memberPCmt(Useed, ren, col, XYsize, circle::Float64, N::Int64)

# Preallocates memory for matrices U and Uaux
Uaux = zeros(ComplexF64, ren, col)
U = zeros(ComplexF64, ren, col)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = Integer.(div.(r.*cos.(theta),XYsize))
yv = Integer.(div.(r.*sin.(theta),XYsize))

# Coherent superposition
@Threads.threads for member=1:N
    Uaux = exp(im * phi[member]) .* MatrixTranslate(Useed, xv[member], yv[member])
    Uaux = Uaux./maximum(abs.(Uaux))
    U = Uaux + U
end

return U
end


function ensemblePC(Useed, Vsize, Hsize, XYsize, circle, N, Ne)

Isum = zeros(Float64, points, points)
Xisum = zeros(Float64, points, points)

    # Generates the desired holograms according to the given parameters
    for memberE = 1:Ne

        # Computes the field of one member
        Upc = memberPC(Useed, Vsize, Hsize, XYsize, circle, N)

        # Incoherent superposition
        I = abs2.(Upc)
        Xi=real(Upc .* conj(rot180(Upc)))
        Isum=Isum+I
        Xisum=Xisum+Xi
    end

    # Averaging
    Iaverage=Isum/Ne
    Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}

end


function ensemblePCmt(Useed, Vsize, Hsize, XYsize, circle, N, Ne)

Isum = zeros(Float64, points, points)
Xisum = zeros(Float64, points, points)

    # Generates the desired holograms according to the given parameters
    for memberE = 1:Ne

        # Computes the field of one member
        Upc = memberPCmt(Useed, Vsize, Hsize, XYsize, circle, N)

        # Incoherent superposition
        I = abs2.(Upc)
        Xi = real(Upc .* conj(rot180(Upc)))
        Isum = Isum+I
        Xisum = Xisum+Xi
    end

    # Averaging
    Iaverage = Isum/Ne
    Xiaverage = Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}

end
