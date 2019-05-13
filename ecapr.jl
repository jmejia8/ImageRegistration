using Metaheuristics

include("tools.jl")
include("appliers.jl")

function myNorm(X, Y)
    sum((Y - X).^2) / size(Y,2)
end

function myError(parameters, x, xDesired)

    y = applyTransformation(x, parameters)

    return myNorm(xDesired, y)

end

function ecapr(X, Y; method=:affine, showResults = true)
    # Y is the desired point set
    η = 2.0
    
    lims= [-1, 1.0]

    if method == :affine
        D = 6
    else
        D = 12
    end

    fitnessFunc(x) = myError(x, X, Y)

    parameters, ee = eca(fitnessFunc, D;η_max = η,
                                        limits = lims,
                                        max_evals=1000D,
                                        correctSol=false, showResults=showResults)

    # performs Y = T(X)
    return parameters, applyTransformation(X, parameters)'
end
