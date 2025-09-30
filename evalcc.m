function [c,flag] = evalcc(n,x,ind)

global problemID pairs nballs a b emin

% ==================================================================

if ( problemID == 1 )

    flag = 1;
   
    c = - x(ind);

    return

end

% ==================================================================

if ( problemID == 2 )

    flag = 1;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    if ( nballs == 1 )
        nck = 0;
    else
        nck = nchoosek(nballs,2);
    end
    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind; 

        c = x.r^2 - ( x.s(i) - 1 )^2 * b^2 * ( cte * x.uv{i}(1)^2 + x.uv{i}(2)^2 );

    elseif ( ind >= nballs + 1 && ind <= nballs + nck )

        i = pairs(ind-nballs, 1);
        j = pairs(ind-nballs, 2);

        c = 4 * x.r^2 - a^2 * ( ( 1 + ( x.s(i) - 1 ) * cte ) * x.uv{i}(1) - ( 1 + ( x.s(j) - 1 ) * cte ) * x.uv{j}(1) )^2 ...
            - b^2 * ( x.s(i) * x.uv{i}(2) - x.s(j) * x.uv{j}(2) )^2;

    elseif ( ind >= nballs + nck + 1 && ind <= 2 * nballs + nck )

        i = ind - ( nballs + nck );

        c = - x.s(i);

    elseif ( ind >= 2 * nballs + nck + 1 && ind <= 3 * nballs + nck )

        i = ind - ( 2 * nballs + nck );

        c = x.s(i) - 1;

    else

        c = - x.r;

    end

    return
end

% ==================================================================


if ( problemID == 3 )

    flag = 1;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    if ind == 1
        A = (x.A + x.A')/2;            % hygiene

        % phi = e(A)^2 = 2*Delta/(s+Delta)  com  p simetrizado
        s     = trace(A);
        s12   = 0.5*(A(1,2)+A(2,1));
        p     = A(1,1)*A(2,2) - s12^2;
        Delta = sqrt(max(s*s - 4*p, 0));
        h     = s + Delta;
        phi   = 2*Delta / max(h,1e-16);

        c = emin^2 - phi;
    else

        flag = 1;
        c = NaN;
    end
    
    return
end