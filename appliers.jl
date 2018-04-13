function applyTransformation(pts, parameters)
    if length(parameters) == 6
        return applyAffine(pts, parameters)
    end

    return applyTransformNoAffine(pts, parameters)
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