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

function applyTransform(x, parameters)

    # rotation and scale matrix
    A = [parameters[1] parameters[2];
         parameters[3] parameters[4]]

    # translate vector
    t = [parameters[5]; parameters[6]]

    cols = size(x, 2)
    return A * x + repmat(t, 1, cols)
end

function applyTransformNoAffine(pts, parameters)
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