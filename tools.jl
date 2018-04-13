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