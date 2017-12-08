using PyPlot
using MATLAB
using Metaheuristics

function plotScatter(x, y, z, w)
    subplot(2,2,1)
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    title("A")
    grid("on")
    
    subplot(2,2,2)
    plot(y[:,1], y[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    title("B")
    grid("on")

    subplot(2,2, 3)
    title("ECA-MP")
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    plot(z[:,1], z[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    grid("on")

    subplot(2,2, 4)
    title("MSAC")
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    plot(w[:,1], w[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    grid("on")

    # xlim(0, 250)
    # ylim(0, 350)
end

function applyAffine(x, parameters)

    # rotation and scale matrix
    A = [parameters[1] parameters[2];
         parameters[3] parameters[4]]

    # translate vector
    t = [parameters[5]; parameters[6]]

    cols = size(x, 2)
    return A * x + repmat(t, 1, cols)
end

function applyQuadratic(pts, parameters)
    p(x,y, a, b, c) = a*x + b*y + c
    g(x,y, a) = a[1]*x.^2 + a[2]*y.^2 + a[3]*x.*y + a[4]*x + a[5]*y + a[6]

    x = pts[1,:]
    y = pts[2,:]

    # a,b,c,d,ee,f = parameters
    # xx = p(x, y, a,  b, c)
    # yy = p(x, y, d, ee, f)
    a = parameters[1:6]
    b = parameters[7:12]

    xx = g(x, y, a)
    yy = g(x, y, b)

    aaa =  [xx yy]'


    return  aaa
end

function applyOrtho(x, parameters)
    # θ: angle
    # T = [x, y]: translate
    # S = [Sx. Sy] scale

    θ = parameters[1]
    T = parameters[2:3]
    S = parameters[4:5]
    # , T, S

    a = cos(-θ)
    b = sin(-θ)
    m = mean(x, 2)

    x = x .- m
    A = [ S[1]*a -b; b S[2]*a ]
    return (A * x) .+ T .+ m


end

function matlabProcrustes(X, Y)
    mat"[d,Z,transform] = procrustes($(X'), $(Y'));
            $c = transform.c;
            $T = transform.T;
            $b = transform.b;
            "
    T *= b
    return  [ reshape(T, 1, 4) c[1,:]' ]
end

function matlabRANSAC(X, Y)
    mat"$Tr = estimateFundamentalMatrix($(X'), $(Y'),...
             'Method','Norm8Point',...
             'NumTrials',20000,'DistanceThreshold',1e-4)"
    
    return [Tr[1], Tr[2],Tr[4], Tr[5],Tr[3], Tr[6]]
end

function matlabAffine(X, Y)
    mat"[tform,inlierPtsDistorted,inlierPtsOriginal] = ...
        estimateGeometricTransform($(Y'),$(X'),...
        'affine');
        $Tr = tform.T
        "

    return [Tr[1], Tr[2],Tr[4], Tr[5],Tr[3], Tr[6]] 
end

function testFuncs(nPoints, i = 1)
    # some trajectories for test algorithm
    circle(θ, r = 1) = r .* cos.(θ), r .* sin.(θ)
    spiral(θ) = exp.(0.1θ).*cos.(4π*θ), exp.(0.1θ).*sin.(4π*θ)
    curve1(θ) = θ.*cos.(θ), sin.(2θ)
    curve2(θ, r) = r.*cos.(2*θ), sin.(2*θ ./ (1+r))


    θ = linspace(0, 2π, nPoints)
    r = linspace(0,  1, nPoints)

    if i == 2
        return spiral(θ)
    elseif i == 3
        return curve1(θ)
    elseif i == 4
        return curve2(θ, r)
    else
        return circle(θ)
    end
    
end

function distanceMatrix(X, Y)
    Nx = size(X, 2)
    Ny = size(Y, 2)

    D = zeros(Ny, Nx)

    for i = 1:Ny
        for j = 1:Nx
            D[i, j] = sum( (X[:,j] - Y[:,i]) .^ 2 )
        end
    end

    return D
end

function hausdorffDistance(X, Y)


    Nx = size(X, 2)
    Ny = size(Y, 2)
    D = distanceMatrix(X, Y)

    h(X, Y, D) = begin
        hh = 0.0
        for i = 1:Ny
            shortest = Inf
            for j = 1:Nx
                if D[i, j] < shortest
                    shortest = D[i, j]
                end
            end

            if shortest > hh
                hh = shortest
            end
        end

        return hh

    end


    return max(h(X, Y, D), h(Y, X, D'))

end