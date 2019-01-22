using Base.Threads   # Support for multithreading
include("/mnt/bendata/Documents/Docs/PhysicsLocal/Beams/Beams.jl")  # include Small-core and LG beams

"""
     PCsimulationSC(N, Ne, w0, wV, circle, TopologicalCharge, points)

Simulate a partially coherent vortex beam using a Small-core expression """
#function PCsimulationSC(N, Ne, w0, wV, circle, TopologicalCharge, points)
function PCsimulationSC(N::Int64, Ne::Int64, w0::Float64, wV::Float64, circle::Float64, TopologicalCharge::Float64, points::Int64)
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

# Incoherent superposition
I = abs2.(U)
Xi=real(U .* conj(rot180(U)))
Isum=Isum+I
Xisum=Xisum+Xi

# Incoherent phase
#println(max(angle.(U)))
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end

"""
     PCsimulationLG(N, Ne, w0, circle, TopologicalCharge, RadialOrder, points)

Simulate a partially coherent vortex beam using a Laguerre-Gaussian expression """
function PCsimulationLG(N::Int64, Ne::Int64; w0::Float64=0.7, circle::Float64=1.0, TopologicalCharge::Int64=1, RadialOrder::Int64=0, points::Int64=128)
# THIS FUNCTION INCORPORATES THE RADIAL ORDER... MAYBE IS SLOWER...
# w0=1
# TopologicalCharge=1
# points=256
# N=10
# Ne=50

# Size of numerical window in terms of circle size
xmax=4*w0

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

for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    # Coherent superposition
#Threads.@threads for member=1:N
    for member=1:N
        Uaux = LaguerreGaussBeam.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], TopologicalCharge, RadialOrder)  # BE CAREFUL WITH THE NORMALIZATION!!!
        Uaux = Uaux/maximum(abs.(Uaux))
        U=Uaux+U
    end

# Incoherent superposition
I = abs2.(U)
Xi=real(U .* conj(rot180(U)))
Isum=Isum+I
Xisum=Xisum+Xi
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end


"""
     PCsimulationHG(N, Ne, w0, circle, TopologicalCharge, RadialOrder, points)

Simulate a partially coherent vortex beam using a Hermite-Gaussian expression """
function PCsimulationHG(N::Int64, Ne::Int64; w0::Float64=0.7, circle::Float64=1.0, m::Int64=1, n::Int64=0, points::Int64=128)
# THIS FUNCTION INCORPORATES THE RADIAL ORDER... MAYBE IS SLOWER...
# w0=1
# TopologicalCharge=1
# points=256
# N=10
# Ne=50

# Size of numerical window in terms of circle size
xmax=3.5*w0

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

for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    # Coherent superposition
Threads.@threads for member=1:N
        Uaux = HermiteGaussBeam.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], m, n)
        Uaux = Uaux/maximum(abs.(Uaux))
        U=Uaux+U
    end

# Incoherent superposition
I = abs2.(U)
Xi=real(U .* conj(rot180(U)))
Isum=Isum+I
Xisum=Xisum+Xi
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end


"""
     PCsimulationIG(N, Ne, w0, circle, p, m, q, kind, points)
p=0,1,2,3... Order of Ince polynomial
0<=m<=p      m is the degree of the Ince polynomial
(p,m) must have the same parity
Simulate a partially coherent vortex beam using a Ince-Gaussian expression """
function PCsimulationIG(N::Int64, Ne::Int64, w0::Float64, circle::Float64, p::Int64, m::Int64, q::Float64, kind::Int64, points::Int64)

# Size of numerical window in terms of circle size
xmax=3.0*w0

# Generates ranges for xs and ys
xs = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
ys = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))

# Generates vector for Matlab-style evaluation
mx, ny = length(xs), length(ys)
Xs = reshape(xs, mx, 1)
Ys = reshape(ys, 1, ny)

# Preallocates memory for matrices Uaux,Isum,Xisum
Uaux = zeros(ComplexF64, points, points)
IGE = zeros(ComplexF64, points, points)
IGO = zeros(ComplexF64, points, points)
Isum = zeros(Float64, points, points)
Xisum = zeros(Float64, points, points)

# Pre-calculates Ince coefficients
if kind == 0
    Cs = CInceCoef(p,m,q)
elseif kind == 1
    Ss = SInceCoef(p,m,q)
else
    Cs = CInceCoef(p,m,q)
    Ss = SInceCoef(p,m,q)
end

for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    # Coherent superposition
    Threads.@threads for member=1:N
        if kind == 0 # IG even
            IGE = IGBeamE.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Cs,))
            #NormN = IGE[abs2.(IGE).==maximum(abs2.(IGE))]
            #NormP = IGE[angle.(IGE).==maximum(angle.(IGE))]
            Uaux = IGE

        elseif kind == 1 # IG odd
            IGO = IGBeamO.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Ss,))
            #NormN = IGO[abs2.(IGO).==maximum(abs2.(IGO))]
            #NormP = IGE[angle.(IGE).==maximum(angle.(IGE))]
            Uaux = IGO

        else # kind == 2 # IG helical
            IGE = IGBeamE.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Cs,))
            NormE = IGE[abs2.(IGE).==maximum(abs2.(IGE))]
            IGE = IGE/NormE[1]
            IGO = IGBeamO.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Ss,))
            NormO = IGO[abs2.(IGO).==maximum(abs2.(IGO))]
            IGO = IGO/NormO[1]
            Uaux = IGE+im*IGO
            #Norm = Uaux[abs2.(Uaux).==maximum(abs2.(Uaux))]
            #Uaux = Uaux/Norm[1]
        end
        U=Uaux+U
    end

# Incoherent superposition
I = abs2.(U)
Xi = real(U .* conj(rot180(U)))
Isum = Isum+I
Xisum = Xisum+Xi
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end


"""
     PCsimulationIGsc(N, Ne, w0, circle, p, m, q, kind, points)
p=0,1,2,3... Order of Ince polynomial
0<=m<=p      m is the degree of the Ince polynomial
(p,m) must have the same parity
Simulate a partially coherent vortex beam using a Ince-Gaussian expression """
function PCsimulationIGsc(N::Int64, Ne::Int64, w0::Float64, circle::Float64, p::Int64, m::Int64, q::Float64, kind::Int64, points::Int64)

# Size of numerical window in terms of circle size
xmax=3.0*w0

# Generates ranges for xs and ys
xs = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
ys = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))

# Generates vector for Matlab-style evaluation
mx, ny = length(xs), length(ys)
Xs = reshape(xs, mx, 1)
Ys = reshape(ys, 1, ny)

# Preallocates memory for matrices Uaux,Isum,Xisum
Uaux = zeros(ComplexF64, points, points)
IGE = zeros(ComplexF64, points, points)
IGO = zeros(ComplexF64, points, points)
Isum = zeros(Float64, points, points)
Xisum = zeros(Float64, points, points)

# Pre-calculates Ince coefficients
if kind == 0
    Cs = CInceCoef(p,m,q)
elseif kind == 1
    Ss = SInceCoef(p,m,q)
else
    Cs = CInceCoef(p,m,q)
    Ss = SInceCoef(p,m,q)
end

for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    for member=1:N
        if kind == 0 # IG even
            IGE = IGBeamE.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Cs,))
            #NormN = IGE[abs2.(IGE).==maximum(abs2.(IGE))]
            #NormP = IGE[angle.(IGE).==maximum(angle.(IGE))]
            Uaux = IGE

        elseif kind == 1 # IG odd
            IGO = IGBeamO.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Ss,))
            #NormN = IGO[abs2.(IGO).==maximum(abs2.(IGO))]
            #NormP = IGE[angle.(IGE).==maximum(angle.(IGE))]
            Uaux = IGO

        else # kind == 2 # IG helical
            IGE = IGBeamE.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Cs,))
            NormE = IGE[abs2.(IGE).==maximum(abs2.(IGE))]
            IGE = IGE/NormE[1]
            IGO = IGBeamO.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], p, m, q, (Ss,))
            NormO = IGO[abs2.(IGO).==maximum(abs2.(IGO))]
            IGO = IGO/NormO[1]
            Uaux = IGE+im*IGO
            #Norm = Uaux[abs2.(Uaux).==maximum(abs2.(Uaux))]
            #Uaux = Uaux/Norm[1]
        end
        U=Uaux+U
    end

# Incoherent superposition
I = abs2.(U)
Xi = real(U .* conj(rot180(U)))
Isum = Isum+I
Xisum = Xisum+Xi
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end


"""
     PCsimulationSUP(N, Ne, w0, circle, TopologicalCharge1, RadialOrder1, TopologicalCharge2, RadialOrder2, points)

Simulate a partially coherent vortex beam using a Laguerre-Gaussian expression """
function PCsimulationSUP(N::Int64, Ne::Int64; w0::Float64=0.7, circle::Float64=1.0, TopologicalCharge1::Int64=1,  RadialOrder1::Int64=0, TopologicalCharge2::Int64=-1, RadialOrder2::Int64=0, points::Int64=128)
# THIS FUNCTION INCORPORATES THE RADIAL ORDER... MAYBE IS SLOWER...
# w0=1
# TopologicalCharge=1
# points=256
# N=10
# Ne=50

# Size of numerical window in terms of circle size
xmax=4*w0

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

Threads.@threads for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    # Coherent superposition
    for member=1:N
        Uaux1 = LaguerreGaussBeam.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], TopologicalCharge1, RadialOrder1)  # BE CAREFUL WITH THE NORMALIZATION!!!
        Uaux1 = Uaux1/maximum(abs.(Uaux1))
        Uaux2 = LaguerreGaussBeam.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], TopologicalCharge2, RadialOrder2)  # BE CAREFUL WITH THE NORMALIZATION!!!
        Uaux2 = Uaux2/maximum(abs.(Uaux2))

        Uaux = Uaux1+Uaux2
        U=Uaux+U
    end

# Incoherent superposition
I = abs2.(U)
Xi=real(U .* conj(rot180(U)))
Isum=Isum+I
Xisum=Xisum+Xi
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end

"""
     PCsimulationMD(N, Ne, w0, circle, TopologicalCharge1, RadialOrder1, TopologicalCharge2, RadialOrder2, points)

Simulate a partially coherent vortex beam using a Laguerre-Gaussian expression """
function PCsimulationMD(N::Int64, Ne::Int64; w0::Float64=0.7, wV::Float64=0.5, circle::Float64=1.0, ltest::Int64=1, TopologicalCharge1::Int64=1,  RadialOrder1::Int64=0, TopologicalCharge2::Int64=-1, RadialOrder2::Int64=0, points::Int64=128)
# THIS FUNCTION INCORPORATES THE RADIAL ORDER... MAYBE IS SLOWER...
# w0=1
# TopologicalCharge=1
# points=256
# N=10
# Ne=50

# Size of numerical window in terms of circle size
xmax=4*w0

# Generates ranges for xs and ys
xs = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
ys = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))

# Generates matrices Matlab-style
mx, ny = length(xs), length(ys)
Xs = reshape(xs, mx, 1)
Ys = reshape(ys, 1, ny)

# Preallocates memory for matrices Uaux,Isum,Xisum
Uaux = zeros(ComplexF64, points, points)
Uaux1 = zeros(ComplexF64, points, points)
Uaux2 = zeros(ComplexF64, points, points)
Umd = zeros(ComplexF64, points, points)
UFFT = zeros(ComplexF64, points, points)
Isum = zeros(Float64, points, points)
Xisum = zeros(Float64, points, points)

for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    # Coherent superposition
    for member=1:N
        Uaux1 = SmallCoreBeam.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], wV, TopologicalCharge1)  # BE CAREFUL WITH THE NORMALIZATION!!!
        Uaux1 = Uaux1/sum(Uaux1)
        #Uaux1 = Uaux1/maximum(abs.(Uaux1))
        Uaux2 = SmallCoreBeam.(Xs .- xv[member], Ys .- yv[member], w0, phi[member], wV, TopologicalCharge2)  # BE CAREFUL WITH THE NORMALIZATION!!!
        Uaux2 = Uaux2/sum(Uaux2)
        #Uaux2 = Uaux2/maximum(abs.(Uaux2))
        Uaux = Uaux1+Uaux2
        Uaux = Uaux/sum(Uaux)  # Normalization
        U=Uaux+U
    end

# Incoherent superposition
#ltest=-1
Umd = exp.(im * ltest * atan2.(Ys,Xs)) .* U
UFFT = fftshift(fft(Umd))

I = abs2.(UFFT)
Xi=real(UFFT .* conj(rot180(UFFT)))
Isum=Isum+I
Xisum=Xisum+Xi
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end



"""
     PCsimulationSUPsc(N, Ne, w0, wV, circle, TopologicalCharge1, TopologicalCharge2, points)

Simulate a partially coherent vortex beam using a Laguerre-Gaussian expression """
function PCsimulationSUPsc(N::Int64, Ne::Int64; w0::Float64=0.7, wV::Float64=1.0, circle::Float64=1.0, TopologicalCharge1::Float64=1, TopologicalCharge2::Float64=-1, points::Int64=128)
# THIS FUNCTION INCORPORATES THE RADIAL ORDER... MAYBE IS SLOWER...
# w0=1
# TopologicalCharge=1
# points=256
# N=10
# Ne=50

# Size of numerical window in terms of circle size
xmax=4*w0

# Generates ranges for xs and ys
xs = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
ys = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))

# Generates matrices Matlab-style
mx, ny = length(xs), length(ys)
Xs = reshape(xs, mx, 1)
Ys = reshape(ys, 1, ny)

# Preallocates memory for matrices Uaux,Isum,Xisum
Uaux = zeros(ComplexF64, points, points)
Isum=zeros(Float64, points, points)
Xisum=zeros(Float64, points, points)

for ensembles=1:Ne  # @threads for multithreading
# Initiliazes to zero
U=zeros(ComplexF64, points, points)

# Uniformly distributed vortices in a circle of radius "circle"
r = circle*sqrt.(rand(N,1))
theta = 2*pi*rand(N,1)
phi = 2*pi*rand(N,1)
xv = r.*cos.(theta)
yv = r.*sin.(theta)

    # Coherent superposition
Threads.@threads for member=1:N
        #Uaux1 = LaguerreGaussBeam.(Xs-xv[member], Ys-yv[member], w0, phi[member], TopologicalCharge1, RadialOrder1)  # BE CAREFUL WITH THE NORMALIZATION!!!
        Uaux1 = SmallCoreBeam.(Xs .- xv[member], Ys .- yv[member], w0, wV, phi[member], TopologicalCharge1)
        Uaux1 = Uaux1/maximum(abs.(Uaux1))
        #Uaux2 = LaguerreGaussBeam.(Xs-xv[member], Ys-yv[member], w0, phi[member], TopologicalCharge2, RadialOrder2)  # BE CAREFUL WITH THE NORMALIZATION!!!
        Uaux2 = SmallCoreBeam.(Xs .- xv[member], Ys .- yv[member], w0, wV, phi[member], TopologicalCharge2)
        Uaux2 = Uaux2/maximum(abs.(Uaux2))

        Uaux = Uaux1+Uaux2
        U=Uaux+U
end

# Incoherent superposition
I = abs2.(U)
Xi=real(U .* conj(rot180(U)))
Isum=Isum+I
Xisum=Xisum+Xi
end

# Averaging
Iaverage=Isum/Ne
Xiaverage=Xisum/Ne

return Iaverage::Array{Float64, 2}, Xiaverage::Array{Float64, 2}
end
