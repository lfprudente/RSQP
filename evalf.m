function [f,flag] = evalf(n,x)

global problemID A Points nPoints LabelPoints

% ==================================================================

if ( problemID == 1 )
    
    flag = 1;
    
    f = - 0.5 * dot(x,A*x);

    return
end

% ==================================================================

if ( problemID == 2 )
    
    flag = 1;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    f = - x.r;

    return
end

% ==================================================================

if ( problemID == 3 )
    
    flag = 1;
    
    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    f = 0;
    for i = 1:nPoints
        fparc = dot(x.A * Points(:,i) + x.b, Points(:,i) ) - 1;

        if ( ( LabelPoints(i) == 1 && fparc > 0 ) || ( LabelPoints(i) == 0 && fparc < 0 ) )
            f = f + fparc^2;
        end
    end
    f =  f / nPoints;

    return
end