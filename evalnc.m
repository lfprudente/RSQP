function [nc,flag] = evalnc(n,x,ind)

global problemID pairs nballs a b N

% ==================================================================

if ( problemID == 1 )

    flag = 1;

    nc = zeros(n,1);

    nc(ind) = - 1;

    return

end

% ==================================================================

if ( problemID == 2 )

    flag = 1;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    nc = zeros(n, 1);

    if ( nballs == 1 )
        nck = 0;
    else
        nck = nchoosek(nballs,2);
    end


    cte = ( b / a )^2;

    if ( ind >= 1 && ind <= nballs )

        i = ind;

        Pui = i;
        Pvi = i + nballs;
        Psi = i + 2 * nballs;
        Pr  = 3 * nballs + 1;

        nc(Pui) = - ( x.s(i) - 1 )^2 * b^2 * cte * 2 * x.uv{i}(1);
        nc(Pvi) = - ( x.s(i) - 1 )^2 * b^2 * 2 * x.uv{i}(2);
        nc(Psi) = - 2 * ( x.s(i) - 1 )* b^2 * ( cte * x.uv{i}(1)^2 + x.uv{i}(2)^2 );
        nc(Pr)  = 2 * x.r;

    elseif ( ind >= nballs + 1 && ind <= nballs + nck )

        i = pairs(ind-nballs, 1);
        j = pairs(ind-nballs, 2);

        Pui = i;
        Pvi = i + nballs;
        Psi = i + 2 * nballs;
        Puj = j;
        Pvj = j + nballs;
        Psj = j + 2 * nballs;
        Pr  = 3 * nballs + 1;

        xxi = ( 1 + ( x.s(i) - 1 ) * cte ) * x.uv{i}(1);
        xxj = ( 1 + ( x.s(j) - 1 ) * cte ) * x.uv{j}(1);
        yyi = x.s(i) * x.uv{i}(2);
        yyj = x.s(j) * x.uv{j}(2);

        nc(Pui) = - 2 * a^2 * ( xxi - xxj ) * ( 1 + ( x.s(i) - 1 ) * cte );
        nc(Pvi) = - 2 * b^2 * ( yyi - yyj ) * x.s(i);
        nc(Psi) = - 2 * a^2 * ( xxi - xxj ) * cte * x.uv{i}(1) - 2 * b^2 * ( yyi - yyj ) * x.uv{i}(2);
        nc(Puj) =   2 * a^2 * ( xxi - xxj ) * ( 1 + ( x.s(j) - 1 ) * cte );
        nc(Pvj) =   2 * b^2 * ( yyi - yyj ) * x.s(j);
        nc(Psj) =   2 * a^2 * ( xxi - xxj ) * cte * x.uv{j}(1) + 2 * b^2 * ( yyi - yyj ) * x.uv{j}(2);
        nc(Pr)  =   8 * x.r;

    elseif ( ind >= nballs + nck + 1 && ind <= 2 * nballs + nck )

        i   = ind - ( nballs + nck );
        Psi = i + 2 * nballs;

        nc(Psi) = - 1;

    elseif ( ind >= 2 * nballs + nck + 1 && ind <= 3 * nballs + nck )

        i   = ind - ( 2 * nballs + nck );
        Psi = i + 2 * nballs;

        nc(Psi) = 1;
    else
        Pr  = 3 * nballs + 1;

        nc(Pr) = - 1;
    end

    return
end

% ==================================================================

if ( problemID == 3 )

    flag = 1;

    if ( ~isstruct(x) )
        [x] = reshapevector(x);
    end

    % Initialize the gradient vector with zeros
    nc = zeros(n, 1);

    if ind == 1
        
        A  = (x.A + x.A')/2;
        s  = trace(A);
        s12= 0.5*(A(1,2)+A(2,1));
        p  = A(1,1)*A(2,2) - s12^2;
    
        Delta = sqrt(max(s*s - 4*p, 0));
        h     = s + Delta;
    
        % phi e first ordem derivatives
        phi   = 2*Delta / max(h,1e-16);
        f     = 1 / max(Delta*h*h, 1e-16);          % f = 1/(Delta*(s+Delta)^2)
    
        phi_s = 8*p * f;
        phi_p = -4*s * f;

        % ds/dA and dp/dA (p parametrized) in order [A11 A21 A12 A22]'
        gs = [1; 0; 0; 1];
        s12 = 0.5*(A(1,2)+A(2,1));
        gp = [A(2,2); -s12; -s12; A(1,1)];
    
        % g5(A) = emin^2 - phi(s,p)  => grad_A g5 = -(phi_s*gs + phi_p*gp)
        gA = -(phi_s*gs + phi_p*gp);
    
        nc(1:4) = gA;                        %
    else
        flag = 1;
        nc(:) = NaN;
    end
end