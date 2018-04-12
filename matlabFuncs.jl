using MATLAB

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