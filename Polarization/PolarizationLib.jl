using Luxor

function PolMap(S0,S1,S2,S3,DenNorm,MapName,sampp=2^4)

    # folder location
    folderloc = "/mnt/bendata/Documents/Docs/PhysicsLocal/SpiralPropagation/"

    # Normalization
    if DenNorm
        deno = S0 .+ 1 .* (S0.==0)
    else
        deno = 1
    end
    s0 = S0./deno
    s1 = S1./deno    # THE MINUS SIGN FIXES A PREVIOUS WRONG SIGN... (stupid frames of reference!)
    s2 = S2./deno
    s3 = S3./deno

    # Geometrical stuff
    A = (s0 .+ hypot.(s1, s2))/2 # Semi-major axis
    B = (s0 .- hypot.(s1, s2))/2 # Semi-minor axis
    PHI = (atan.(s2,s1))./2       # Angle
    HD = sign.(s3);           # Handedness

    # Size of numerical window in terms of circle size
    w0 = 0.5
    xmax = 2*w0
    points = 512

    # Generates ranges for xs and ys
    xs = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
    ys = xmax*(2/points)*collect(range(-points/2,length=points,stop=points/2-1))
    Xs = reshape(xs,points,1)
    Ys = reshape(ys,1,points)

    # Samples from the matrices that define the ellipses
    space = div(points,sampp)
    As = A[space+1:space:end, space+1:space:end]
    Bs = B[space+1:space:end, space+1:space:end]
    PHIs = PHI[space+1:space:end, space+1:space:end]
    HDs = HD[space+1:space:end, space+1:space:end]
    ECCs = sqrt.(1 .- (Bs.^2) ./ (As.^2))

    # get coordinates
    Xsel = repeat(Xs[space+1:space:end],1,sampp-1)
    Ysel = repeat(Ys[:,space+1:space:end],sampp-1,1)

    # Scaling parameters
    scalingSize = 12 + 10*(DenNorm==false) - (0.1 - 0.5*(DenNorm==false))*sampp
    scalingSpace = 10e1

    Drawing(2*scalingSpace*xmax,2*scalingSpace*xmax,
            "$(folderloc)$(MapName)")
    sethue("red")
    origin()

    # "Fixes" coordinate system
    transform([-1 0 0 1 0 0])
    transform([cos(pi), -sin(pi), sin(pi), cos(pi), 0, 0])
    #rulers()

    # Plots ellipses
    for ii = 1:sampp-1
        for jj = 1:sampp-1
            gsave()
            setline(0.05)
            Luxor.translate(scalingSpace*Xsel[ii, jj], scalingSpace*Ysel[ii, jj])
            Luxor.rotate(PHIs[ii,jj])
            if HDs[ii,jj] > 0
                sethue("blue")
            else
                sethue("red")
            end
            #sethue("black")
            setline(log(0.2+ECCs[ii,jj])+0.1)   # Scaling width of the line
            ellipse(0, 0, scalingSize*As[ii, jj], scalingSize*Bs[ii, jj], :fillstroke)
            grestore()
        end
    end
    finish()
    preview()
end

function Halfwaveplate(Ex, Ey, theta)
    Jhwp = [cos(2*theta) sin(2*theta); sin(2*theta) -cos(2*theta)]
    EX = Jhwp[1,1]*Ex + Jhwp[1,2]*Ey
    EY = Jhwp[2,1]*Ex + Jhwp[2,2]*Ey
    return EX, EY
end

function Quarterwaveplate(Ex, Ey, theta)
    Jqwp = (1/sqrt(2))*[1+im*cos(2*theta) im*sin(2*theta); im*sin(2*theta) 1-im*cos(2*theta)];
    EX = Jqwp[1,1]*Ex .+ Jqwp[1,2]*Ey
    EY = Jqwp[2,1]*Ex .+ Jqwp[2,2]*Ey
    return EX, EY
end
