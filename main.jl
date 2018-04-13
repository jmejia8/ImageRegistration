using Metaheuristics

include("tools.jl")
include("appliers.jl")
include("matlabFuncs.jl")
include("createFigures.jl")

function myNorm(X, Y)
    sum((Y - X).^2) / size(Y,2)
end

function myError(parameters, x, xReal)

    y = applyTransformation(x, parameters)

    return myNorm(xReal, y)

end

function main()
    Y = readcsv("data/test6.csv")'
    X = readcsv("data/test6-1.csv")'

    η = 2.0
    
    lims= [-5, 5.0]
    D = 6

    fitnessFunc(x) = myError(x, Y, X)

    @time parameters, ee = eca(fitnessFunc, D;η_max = η,
                                        limits = lims,
                                        max_evals=1000D,
                                        # termination= x->std(-1 + 1.0./x) < 1e-10,
                                        correctSol=false)

    X_approx = applyTransformation(Y, parameters)
    
    @time Tr = matlabAffine(X, Y)
    X_matlab = applyTransformation(Y, Tr)
    plotScatter(X', Y', X_approx', X_matlab')

    @printf("ECA  error: %e\n", myNorm(X, X_approx))
    @printf("MLAB error: %e\n", myNorm(X, X_matlab))

    println(parameters)
end

main()