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

function myNorm(X, Y)
    sum((Y - X).^2) / size(Y,2)
end

function myError(parameters, x, xReal)

    y = applyTransformNoAffine(x, parameters)
    # err = norm(y - xReal)
    err = myNorm(xReal, y)# sum((y - xReal).^2) / size(y,2)


    return err


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

function testRW()
    Y = readcsv("data/test6.csv")'
    X = readcsv("data/test6-1.csv")'

    η = 5.0
    
    lims= (-100, 100)
    D = 6

    fitnessFunc(x) = myError(x, Y, X)

    @time parameters, ee = eca(fitnessFunc, D;η_max = η,
                                        limits = lims,
                                        max_evals=100000D,
                                        termination= x->std(-1 + 1.0./x) < 1e-10,
                                        correctSol=false)

    X_approx = applyTransformNoAffine(Y, parameters)
    
    @time Tr = matlabAffine(X, Y)
    X_matlab = applyTransform(Y, Tr)
    plotScatter(X', Y', X_approx', X_matlab')

    @printf("ECA  error: %e\n", myNorm(X, X_approx))
    @printf("MLAB error: %e\n", myNorm(X, X_matlab))

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

function main()
    # number of points
    nPoints = 100

    x, y = testFuncs(nPoints, 1)

    # Data
    X = 10 + [x y]'


    # parameters for ECA algorithm
    η   = 5.0
    lims= (-100, 100)
    D   = 12

    for i in 1:1
        # Random affine transformation
        originalParms = -100 + 200rand(12)
        fitnessFunc(x) = myError(x, Y, X)

        # uncomment for additive noise
        Y = applyTransformNoAffine(X, originalParms )
        Y = applyTransformNoAffine(Y, originalParms )

        # Find affine transformation
        @time approxParms, ee = eca(fitnessFunc, D; η_max = η,
                                                    limits = lims,
                                                    max_evals=100000D,
                                                    termination= x->std(-1 + 1.0./x) < 1e-10,
                                                    correctSol=false)

        X_approx = applyTransformNoAffine(Y, approxParms)

        # Calculate affine transformation using MATLAB
        @time Tr = matlabAffine(X, Y)
        X_matlab = applyTransform(Y, Tr)

        # Plot results
        plotScatter(X', Y', X_approx', X_matlab')

        # comparing error
        @printf("ECA  error: %e\n", myNorm(X, X_approx))
        @printf("MLAB error: %e\n", myNorm(X, X_matlab))
    end
end

main()
# testRW()