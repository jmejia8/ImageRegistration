function test()
    # number of points
    nPoints = 100

    x, y = testFuncs(nPoints, 1)

    # Data
    X = 10 + [x y]'


    # parameters for ECA algorithm
    η   = 5.0
    lims= [-100.0, 100]
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
                                                    max_evals=10000D,
                                                    # termination= x->std(-1 + 1.0./x) < 1e-10,
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