using Base.Threads   # Support for multithreading
using FFTW
include("Beams.jl")  # include Small-core and LG beams

"""
     PCsimulationPhase(N, Ne, w0, wV, circle, TopologicalCharge, points)

Computes the "phase" of a partially coherent field """
function PCsimulationPhase(N::Int64, Ne::Int64, w0::Float64, wV::Float64, circle::Float64, TopologicalCharge::Float64, points::Int64)
# w0=1
# TopologicalCharge=1
# points=256
# N=10
# Ne=50

# Size of numerical window
xmax=2.5*w0

# Generates ranges for xs and ys
xs = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
ys = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))

# Generates matrices Matlab-style
mx, ny = length(xs), length(ys)
Xs = reshape(xs, mx, 1)
Ys = reshape(ys, 1, ny)

# Preallocates memory for matrices Uaux,Isum,Xisum
Uaux = zeros(ComplexF64, points, points)
Ph = zeros(Float64, points, points)
Phsum = zeros(Float64, points, points)
Xisum = zeros(Float64, points, points)


#Threads.@threads
for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
#U= zeros(Complex, points, points)
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    # Coherent superposition
    Threads.@threads for member=1:N
            Uaux = SmallCoreBeam.(Xs .- xv[member], Ys .- yv[member], w0, wV, phi[member], TopologicalCharge)
            Uaux = Uaux/maximum(abs.(Uaux))
            U=Uaux+U
    end

# Incoherent superposition
I = abs2.(U)
#Xi=real(U .* conj(rot180(U)))
Ph = angle.(U)
#PhXi = angle.(Xi)
#Isum=Isum+I

# "Incoherent" phase
Phsum=Phsum+Ph
end

# Averaging
Phaverage=Phsum/Ne

return Ph::Array{Float64, 2}, Phsum::Array{Float64, 2}, Phaverage::Array{Float64, 2}
end


"""
     ModalDecomp(N, Ne, w0, wV, circle, TopologicalCharge, points)

Modal decomposition with PCVB """
function ModalDecomp(N::Int64, Ne::Int64, w0::Float64, wV::Float64, circle::Float64, TopologicalCharge::Float64, points::Int64, PhaseModal::Array{Float64, 2})

# numerical window
xmax=2.5*w0

# Generates ranges for xs and ys
xs = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
ys = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))

# Generates matrices Matlab-style
mx, ny = length(xs), length(ys)
Xs = reshape(xs, mx, 1)
Ys = reshape(ys, 1, ny)

# Preallocates memory for matrices Uaux,Isum,Xisum
Uaux = zeros(ComplexF64, points, points)
Isum = zeros(Float64, points, points)
Xisum = zeros(Float64, points, points)

#Threads.@threads
for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
#U= zeros(Complex, points, points)
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

# Coherent superposition
Threads.@threads for member=1:N
        Uaux = SmallCoreBeam.(Xs .- xv[member], Ys .- yv[member], w0, wV, phi[member], TopologicalCharge)
        Uaux = Uaux/maximum(abs.(Uaux))
        U=Uaux+U
end

Umd = U .* exp.(im * PhaseModal)

# Phase element for modal de
UFFT = fftshift(fft(Umd))

# Incoherent superposition
I = abs2.(UFFT)
#Xi=real(U .* conj(rot180(U)))

Isum=Isum+I
#Xisum=Xisum+Xi

# Incoherent phase
#println(max(angle.(U)))
end

# Averaging
Iaverage=Isum/Ne
#Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}
end
