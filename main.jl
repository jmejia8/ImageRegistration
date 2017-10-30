using PyPlot
using MATLAB
using Metaheuristics

function plotScatter(x, y, z, w)
    subplot(1,3,1)
    title("Transformation")
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    plot(y[:,1], y[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    grid("on")

    subplot(1,3, 2)
    title("ECA")
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    plot(z[:,1], z[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    grid("on")

    subplot(1,3, 3)
    title("MATLAB")
    plot(x[:,1], x[:,2], marker=:o, lw = 0, markersize=4, color=:b)
    plot(w[:,1], w[:,2], marker=:o, lw = 0, markersize=4, color=:r)
    grid("on")

    # xlim(0, 250)
    # ylim(0, 350)
end

function applyTransform(x, parameters)
    # A = [prms1 prms2; prms3 prms4;]
    # b = [prms4, prms6]

    # rotation and scale matrix
    A = [parameters[1] parameters[2];
         parameters[3] parameters[4]]

    # translate vector
    t = [parameters[5]; parameters[6]]

    cols = size(x, 2)
    return A * x + repmat(t, 1, cols)
end

function a_pplyTransform(x, parameters)
    p(t, a, b, c) = a* exp.(b * t) + c  


    a,b,c,d,ee,f = parameters
    xx = p(x[1,:], a,  b, c)
    yy = p(x[2,:], d, ee, f)
    aaa =  [xx  yy]'


    return  aaa
end

function myError(population, x, xReal)
    rows, cols = size(population, 1, 2)

    if cols == 1
        population = population'
        rows, cols = size(population, 1, 2)
    end

    errors = zeros(rows)
    for i = 1:rows
        y = applyTransform(x, population[i,:])
        errors[i] = norm(y - xReal)
    end

    if rows == 1
        return errors[1]
    end

    return errors


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
    # reshape(Tr[1:2,:], 1, 6)
end

function matlabAffine(X, Y)
    mat"[tform,inlierPtsDistorted,inlierPtsOriginal] = ...
        estimateGeometricTransform($(Y'),$(X'),...
        'affine');
        $Tr = tform.T
        "

    return [Tr[1], Tr[2],Tr[4], Tr[5],Tr[3], Tr[6]] 
    # reshape(Tr[1:2,:], 1, 6)
end

function testRW()
    X = readcsv("data/test1.csv")'
    Y = readcsv("data/test1-2.csv")'

    K     = 7
    η_max = 5.0
    
    limits= (-100, 100)
    D = 6
    N = K*D

    max_evals = 10000D
    fitnessFunc(x) = myError(x, Y, X)

    parameters, ee = eca(fitnessFunc, D, N, max_evals, K, η_max, limits, 0)

    X_approx = applyTransform(Y, parameters)
    Tr = matlabAffine(X, Y)
    X_matlab = applyTransform(Y, Tr)
    plotScatter(X', Y', X_approx', X_matlab')

    @printf("ECA  error: %e\n", norm(X - X_approx))
    @printf("MLAB error: %e\n", norm(X - X_matlab))

end


function main()
    circle(θ, r = 1) = r .* cos.(θ), r .* sin.(θ)
    spiral(θ) = exp.(0.1θ).*cos.(4π*θ), exp.(0.1θ).*sin.(4π*θ)
    curve1(θ) = θ.*cos.(θ), sin.(2θ)
    curve2(θ, r) = r.*cos.(2*θ), sin.(2*θ ./ (1+r))

    nPoints = 100

    θ = linspace(0, 2π, nPoints)
    r = linspace(0,  1, nPoints)
    
    x, y = curve1(θ)

    X = 10 + [x y]'


    η = 5.0
    
    lims= (-1, 1)
    D = 10

    max_evals = 10000D

    for i in 1:5
        originalParms = 10randn(6)
        fitnessFunc(x) = myError(x, Y, X)

        Y = applyTransform(X, originalParms)  #+ randn(2, nPoints)


        @time approxParms, ee = eca(fitnessFunc, D; η_max = η,
                                                    limits = lims,
                                                    termination= x->std(x) < 1e-15,
                                                    correctSol=false)

        X_approx = applyTransform(Y, approxParms)

        @time Tr = matlabAffine(X, Y)
        X_matlab = applyTransform(Y, Tr)
        plotScatter(X', Y', X_approx', X_matlab')

        @printf("ECA  error: %e\n", norm(X - X_approx))
        @printf("MLAB error: %e\n", norm(X - X_matlab))
    end
end

main()
# testRW()