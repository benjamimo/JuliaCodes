"""
     HoloCream(U, lines, alf)

Creates the hologram of the desired field U """
function HoloCrea(U, lines, alf)

    # SPATIAL SCALE
    Ux=15.36E-3;
    Uy=8.64E-3;
    Hsize=size(U,2);
    Vsize=size(U,1);

    # Generates ranges for xs and ys
    ys = Ux*(2/Hsize)*collect(range(-Hsize/2,length=Hsize,stop=Hsize/2-1))
    xs = Uy*(2/Vsize)*collect(range(-Vsize/2,length=Vsize,stop=Vsize/2-1))

    # Coordinates
    mx, ny = length(xs), length(ys)
    Ys = reshape(xs, mx, 1)
    Xs = reshape(ys, 1, ny)

    # Hologram creation: receives complex field U properly normalized such that U is in [0,1]
    A = abs.(U)
    g = A./maximum(A)
    phase = angle.(U)

    # lines = 100
    linespMM = 2*lines/(1E-2)
    # alf = 0
    kxs = linespMM*cos(alf)
    kys = linespMM*sin(alf)
    H = g .*(mod.(phase .+ (kxs*Xs .+ kys*Ys), 2*pi)./pi .- 1) .+ 1;
    H = H./maximum(H);

    return H
end
