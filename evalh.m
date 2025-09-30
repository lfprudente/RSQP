function [h,flag] = evalh(n,x)

global problemID A Points nPoints LabelPoints

% ==================================================================

if ( problemID == 1 )
    
    flag = 1;
    
    h = - A;

    return
end

% ==================================================================

if ( problemID == 2 )
    
    flag = 0;
    
    h = zeros(n);

    return
end

% ==================================================================

if ( problemID == 3 )

    flag = 1;

    if ~isstruct(x)
        x = reshapevector(x); % user-provided reshaper
    end

    h = zeros(6,6);

    % -- Build Hessian as sum_i 2 * c_i c_i' / nPoints over active i --
    for i = 1:nPoints
        pi1 = Points(1,i);
        pi2 = Points(2,i);

        % f_i(x) = (A p_i + b)' p_i - 1
        f_i = (x.A(1,1)*pi1 + x.A(1,2)*pi2 + x.b(1))*pi1 + ...
              (x.A(2,1)*pi1 + x.A(2,2)*pi2 + x.b(2))*pi2 - 1;

        % Active set exactly as in evalf/evalg:
        is_active = (LabelPoints(i) == 1 && f_i > 0) || (LabelPoints(i) == 0 && f_i < 0);
        if ~is_active, continue; end

        % grad f_i wrt [A11 A21 A12 A22 b1 b2]
        c = [pi1^2; pi1*pi2; pi1*pi2; pi2^2; pi1; pi2];  % 6x1

        % Hessian contribution: 2 * c c^T
        h = h + 2*(c*c.');
    end

    h = h / nPoints;
end