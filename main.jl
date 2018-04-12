include("tools.jl")



function myNorm(X, Y)
    return sum( (Y - X).^2 ) / size(Y,2)
end

# function myNorm(X, Y)
#     Nx = size(X, 2)
#     Ny = size(Y, 2)

#     μ = zeros(Int64, Ny, Nx)
#     s = 0.0
#     for i = 1:Ny

#         m = -1.0
#         jmin = 0
#         for j= 1:Nx
#             if sum(μ[:,j]) != 0
#                 continue
#             end

#             d = norm(X[:, j] - Y[:, i])

#             if m < 0 || d < m
#                 m = d
#                 jmin = j
#             end

#         end
#         μ[i, jmin] = 1
#         s += m
#     end

#     return  s / Ny
# end

function myError(parameters, X, Y, method, T)
    # if method == "affine"
    #     return myNorm(X, applyAffine(Y, parameters))
    # elseif method == "ortho"
    #     return myNorm(X, applyOrtho(Y, parameters))
    # elseif method == "quadratic"
    #     return myNorm(X, applyQuadratic(Y, parameters))
    # end

    myNorm(X, T(Y, parameters))


end

function strategy(method = "affine")
    # methods: ortho, affine, quadratic
    η = 5.0
    lims = [-2, 2]
    D = 6
    correctSolution = false
    func = applyAffine

    if method == "ortho"
        η = 2.0
    
        lims= [ 0.0 -100.0  -10.0 -5.0 -5.0; 
                2π   100     10   5.0 5.0]
        D = 5
        correctSolution = true
        func = applyOrtho
    elseif method == "quadratic"
        lims = [-2, 2]
        D = 12
        func = applyQuadratic
    elseif method == "quadratic2"
        lims = [-2, 2]
        func = applyQuadratic2
    end

    return η, lims, D, correctSolution, func
end


function testRW()
    nmes = ["./convergence_sin_affine.csv",
            "./convergence_chin_affine.csv",
            "./convergence_papa_affine.csv",
            "./convergence_papa2_affine.csv",
            "./convergence_pez_affine.csv",
           ]

zz = 5
for i = 1:zz

    X = readcsv("data/test$(i+1).csv")'
    Y = readcsv("data/test$(i+1)-1.csv")'

    method = "quadratic"
    η, lims, D, correctSolution, T = strategy(method)

    fitnessFunc(x, method = method, Tr = T) = myError(x, X, Y, method, Tr)

    @time parameters, ee = eca(fitnessFunc, D;η_max = η,
                                        N = 10D,
                                        limits = lims,
                                        max_evals=10000D,
                                        termination= x->std(-1.0 + 1.0 ./ x) < 1e-10,
                                        saveConvergence = true,
                                        correctSol=correctSolution)

    println(parameters)

    X_approx = T(Y, parameters)
    
    @time Tr = matlabAffine(X, Y)
    X_matlab = applyAffine(Y, Tr)
    @printf("MSAC error: %e\n", myNorm(X, X_matlab))
    continue
    # plotScatter(X', Y', X_approx', X_matlab')

    subplot(5, 3, 3i-2)
    plot(X[1,:], X[2,:], marker=:o, lw=0, color="blue", markersize=4)
    plot(Y[1,:], Y[2,:], marker=:o, lw=0, color="red", markersize=4)
    grid("on")

    subplot(5, 3, 3i-1)
    plot(X[1,:], X[2,:], marker=:o, lw=0, color="blue", markersize=4)
    plot(X_approx[1,:], X_approx[2,:], marker=:o, lw=0, color="red", markersize=4)
    grid("on")

    subplot(5, 3, 3i)
    a = readcsv(nmes[i])

    b = 1:length(a)
    a = 10 + log.(1.0 - a)

    plot(b, a, color="black")
    xlabel("Generation")
    ylabel("Log Error")
    title("Convergence")
    grid("on")

    # @printf("ECA  error: %e\n", myNorm(X, X_approx))
    @printf("MSAC error: %e\n", myNorm(X, X_matlab))
end
end



function main()
    # number of points
    nPoints = 50

    x, y = testFuncs(nPoints, 5)

    # Data
    X = 10 + [x y]'


    # parameters for ECA algorithm
    method = "affine"
    η, lims, D, correctSolution, T = strategy(method)

    for i in 1:1
        # Random affine transformation
        fitnessFunc(x, method = method, Tr = T) = myError(x, X, Y, method, Tr)

        # uncomment for additive noise
        originalParms = 10randn(12)
        Y = applyQuadratic(X, originalParms )
        Y = applyQuadratic(Y, -originalParms )

        # writecsv("points.csv", Y)

        # Find affine transformation
        @time approxParms, ee = eca(fitnessFunc, D; η_max = η,
                                                        N = 10D,
                                                    limits = lims,
                                                    max_evals=30000D,
                                                    saveGens = false,
                                                    termination= x->std( -1.0 + 1.0 ./ x) < 1e-10,
                                                    correctSol=false)

        X_approx = T(Y, approxParms)

        # Calculate affine transformation using MATLAB
        @time Tr = matlabAffine(X, Y)
        X_matlab = applyAffine(Y, Tr)

        # Plot results
        plotScatter(X', Y', X_approx', X_matlab')

        # comparing error
        @printf("ECA  error: %e\n", myNorm(X, X_approx))
        @printf("MSAC error: %e\n", myNorm(X, X_matlab))
    end
end

# main()
testRW()