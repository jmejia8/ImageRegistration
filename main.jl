using Metaheuristics

include("tools.jl")
include("appliers.jl")
include("matlabFuncs.jl")
include("createFigures.jl")

function myNorm(X, Y)
    sum((Y - X).^2) / size(Y,2)
end

function myError(parameters, x, xDesired)

    y = applyTransformation(x, parameters)

    return myNorm(xDesired, y)

end

function ecapr(X, Y, method=:affine)
    # Y is the desired point set
    η = 2.0
    
    lims= [-5, 5.0]

    if method == :affine
        D = 6
    else
        D = 12
    end

    fitnessFunc(x) = myError(x, X, Y)

    parameters, ee = eca(fitnessFunc, D;η_max = η,
                                        limits = lims,
                                        max_evals=1000D,
                                        correctSol=false)

    # performs Y = T(X)
    return applyTransformation(X, parameters)
end

function main()
    # observed
    X = readcsv("data/test6.csv")'

    # desired
    Y = readcsv("data/test6-1.csv")'

    Y_ecapr = ecapr(X, Y)

    Y_msac = applyTransformation(X, matlabAffine(X, Y))

    plotScatter(Y', X', Y_ecapr', Y_msac')

    @printf("ECA  error: %e\n", myNorm(Y, Y_ecapr))
    @printf("MSAC error: %e\n", myNorm(Y, Y_msac))

end

main()