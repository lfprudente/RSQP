function [g,flag] = evalg(n,x)

global problemID A Points nPoints LabelPoints

% ==================================================================

if ( problemID == 1 )
    
    flag = 1;

    g = zeros(n,1);

    g = - A * x;

    return
end

% ==================================================================

if ( problemID == 2 )
    
    flag = 1;
    
    g = zeros(n, 1);

    g(end) = - 1;

    return
end

% ==================================================================

if ( problemID == 3 )
    
    flag = 1;
    
    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    g = zeros(n,1);

    for i = 1:nPoints
        fparc = dot(x.A * Points(:,i) + x.b, Points(:,i) ) - 1;

        if ( ( LabelPoints(i) == 1 && fparc > 0 ) || ( LabelPoints(i) == 0 && fparc < 0 ) )
            g(1) = g(1) + 2 * fparc * Points(1,i)^2;
            g(2) = g(2) + 2 * fparc * Points(1,i) *  Points(2,i);
            g(3) = g(3) + 2 * fparc * Points(1,i) *  Points(2,i);
            g(4) = g(4) + 2 * fparc * Points(2,i)^2;
            g(5) = g(5) + 2 * fparc * Points(1,i);
            g(6) = g(6) + 2 * fparc * Points(2,i);
        end
    end
    
    g =  g / nPoints;

    return
end